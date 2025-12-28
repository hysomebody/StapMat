% TETRAELEMENT 4-Node Tetrahedral Element Class (Linear)
%
% 功能：实现四节点四面体单元
%   适配新架构：每个节点 6 个自由度，但在刚度矩阵中仅填充平动 (1,2,3) 部分。
%
% Call procedures:
%   ./Element.m
%   ./Node.m
%   ./Material.m (Expects properties E, Nu, Rho)
%
%
% By 程云志

classdef TetraElement < Element
    methods
        
        function obj = TetraElement()
            obj@Element();
        end

       
        function Read(obj, fid, expectedID, materialList, globalNodeList)
            lineStr = fgetl(fid);
            data = str2num(lineStr); %#ok<ST2NM>

            
            inputID = round(data(1));
            if inputID ~= expectedID
                error('Element ID mismatch. Expected: %d, Read: %d', expectedID, inputID);
            end
            obj.ID = inputID;

            
            nIndices = [round(data(2)), round(data(3)), round(data(4)), round(data(5))];
            obj.SetNodes(globalNodeList, nIndices);

            
            matID = round(data(6));
            obj.Material = materialList(matID);
        end

        % =================================================================
        % CalcStiffness (计算刚度矩阵)
        % Output: k (12x12 matrix, 4 nodes * 3 DOFs)
        % =================================================================
        function k = CalcStiffness(obj)
            % 1. Get Coordinates
            x = zeros(4,1); y = zeros(4,1); z = zeros(4,1);
            for i = 1:4
                node = obj.Nodes(i);
                x(i) = node.XYZ(1);
                y(i) = node.XYZ(2);
                z(i) = node.XYZ(3);
            end

            % 2. Calculate Volume (V6 = 6 * Volume)
            M_vol = [1, x(1), y(1), z(1);
                     1, x(2), y(2), z(2);
                     1, x(3), y(3), z(3);
                     1, x(4), y(4), z(4)];
            V6 = det(M_vol);
            V = abs(V6) / 6.0;

            % 3. Calculate Strain-Displacement Matrix B (6x12)
            B = zeros(6, 12);
            
            % i: Node index (1-4)
            for i = 1:4
                j = mod(i, 4) + 1;
                m = mod(i + 1, 4) + 1;
                p = mod(i + 2, 4) + 1;

                % Cofactors
                bi = -det([1, y(j), z(j); 1, y(m), z(m); 1, y(p), z(p)]);
                ci =  det([1, x(j), z(j); 1, x(m), z(m); 1, x(p), z(p)]);
                di = -det([1, x(j), y(j); 1, x(m), y(m); 1, x(p), y(p)]);

                % Fill B matrix blocks (3 columns per node)
                col_start = (i-1)*3 + 1;
                col_end   = (i-1)*3 + 3;
                
                % Standard FEM B-matrix structure for 3D Solid
                % B_i = [bi 0 0; 0 ci 0; 0 0 di; 0 di ci; di 0 bi; ci bi 0]
                B(:, col_start:col_end) = (-1)^(i+1) * ...
                    [bi,  0,  0;
                      0, ci,  0;
                      0,  0, di;
                      0, di, ci;
                     di,  0, bi;
                     ci, bi,  0];
            end
            B = B / V6;

            % 4. Get Constitutive Matrix D (6x6)
            D = obj.GetDMatrix();

            % 5. Compute Element Stiffness (12x12)
            % k = integral(B' * D * B) dV = B' * D * B * V
            k = B' * D * B * V;
            
        end

        function f = CalcThermalLoad(obj)
            % 1. 获取节点坐标
            x=zeros(4,1); y=zeros(4,1); z=zeros(4,1);
            temps = zeros(4,1);
            for i=1:4 
                n=obj.Nodes(i); 
                x(i)=n.XYZ(1); y(i)=n.XYZ(2); z(i)=n.XYZ(3);
                % 假设 Node 类有 Temperature 属性
                if isprop(n, 'Temperature')
                    temps(i) = n.Temperature;
                else
                    temps(i) = 0;
                end
            end

            % 2. 计算体积 V6 和 V
            M_vol = [1, x(1), y(1), z(1); 
                     1, x(2), y(2), z(2); 
                     1, x(3), y(3), z(3); 
                     1, x(4), y(4), z(4)];
            V6 = det(M_vol);
            V = abs(V6) / 6.0;

            % 3. 计算 B 矩阵 (6x12)
            B = zeros(6, 12);
            for i = 1:4
                j=mod(i,4)+1; m=mod(i+1,4)+1; p=mod(i+2,4)+1;
                bi = -det([1, y(j), z(j); 1, y(m), z(m); 1, y(p), z(p)]);
                ci =  det([1, x(j), z(j); 1, x(m), z(m); 1, x(p), z(p)]);
                di = -det([1, x(j), y(j); 1, x(m), y(m); 1, x(p), y(p)]);
                
                col_start = (i-1)*3 + 1;
                col_end   = (i-1)*3 + 3;
                
                B(:, col_start:col_end) = (-1)^(i+1) * ...
                    [bi,  0,  0;
                      0, ci,  0;
                      0,  0, di;
                      0, di, ci;
                     di,  0, bi;
                     ci, bi,  0];
            end
            B = B / V6;

            % 4. 获取 D 矩阵
            D = obj.GetDMatrix();
            
            % 5. 计算热载荷
            % 假设 Alpha 存在于 Material 中
            alpha = 0;
            if isprop(obj.Material, 'Alpha')
                alpha = obj.Material.Alpha;
            end
            
            % 计算平均温升 (假设初始温度为0，即 dT = T_current)
            T_avg = sum(temps) / 4.0;
            
            % 热应变向量 [alpha*dT, alpha*dT, alpha*dT, 0, 0, 0]'
            eps_th = alpha * T_avg * [1; 1; 1; 0; 0; 0];
            
            % 计算等效节点力: f = Integral(B' * D * eps_th) dV
            % 对于常应变单元，积分为直接乘以体积 V
            f = V * (B' * D * eps_th);
        end

        % =================================================================
        % CalcStress (计算单元应力 - 修正版含热应力)
        % =================================================================
        function results = CalcStress(obj, globalU)
            % 1. Extract Local Displacements (Translational Only)
            u_elem = zeros(12, 1);
            for i = 1:4
                node = obj.Nodes(i);
                % Get U, V, W from global vector
                for dof = 1:3
                    eqID = node.BCode(dof);
                    if eqID > 0
                        u_elem((i-1)*3 + dof) = globalU(eqID);
                    end
                end
            end

            % 2. Re-calculate B and V6 (Needed for strain)
            x = zeros(4,1); y = zeros(4,1); z = zeros(4,1);
            temps = zeros(4,1);
            for i = 1:4
                n = obj.Nodes(i); 
                x(i)=n.XYZ(1); y(i)=n.XYZ(2); z(i)=n.XYZ(3);
                if isprop(n, 'Temperature'), temps(i) = n.Temperature; end
            end
            
            M_vol = [1, x(1), y(1), z(1); 1, x(2), y(2), z(2); 1, x(3), y(3), z(3); 1, x(4), y(4), z(4)];
            V6 = det(M_vol);
            
            B = zeros(6, 12);
            for i = 1:4
                j=mod(i,4)+1; m=mod(i+1,4)+1; p=mod(i+2,4)+1;
                bi = -det([1, y(j), z(j); 1, y(m), z(m); 1, y(p), z(p)]);
                ci =  det([1, x(j), z(j); 1, x(m), z(m); 1, x(p), z(p)]);
                di = -det([1, x(j), y(j); 1, x(m), y(m); 1, x(p), y(p)]);
                B(:,(i-1)*3+1:(i-1)*3+3) = (-1)^(i+1)*[bi,0,0; 0,ci,0; 0,0,di; 0,di,ci; di,0,bi; ci,bi,0];
            end
            B = B / V6;

            % 3. Calculate Total Strain
            % Strain = B * u
            strain_total = B * u_elem;
            
            % 4. Calculate Thermal Strain
            % 获取热膨胀系数
            alpha = 0;
            if isprop(obj.Material, 'Alpha')
                alpha = obj.Material.Alpha;
            end
            
            % 平均温度
            T_avg = sum(temps) / 4.0;
            
            % 热应变向量 [alpha*dT, alpha*dT, alpha*dT, 0, 0, 0]'
            eps_th = alpha * T_avg * [1; 1; 1; 0; 0; 0];

            % 5. Calculate Mechanical Stress
            D = obj.GetDMatrix();
            % Stress = D * (Strain_Total - Strain_Thermal)
            % [关键修改]：这里必须减去 eps_th
            stress_vec = D * (strain_total - eps_th); 
            
            % stress_vec format: [sig_x, sig_y, sig_z, tau_yz, tau_zx, tau_xy]'

            % 6. Calculate Von Mises Stress
            sx = stress_vec(1); sy = stress_vec(2); sz = stress_vec(3);
            t_yz = stress_vec(4); t_zx = stress_vec(5); t_xy = stress_vec(6);

            % Von Mises Formula:
            vm_sq = 0.5 * ((sx-sy)^2 + (sy-sz)^2 + (sz-sx)^2 + 6*(t_xy^2 + t_yz^2 + t_zx^2));
            von_mises = sqrt(vm_sq);

            results.StressVector = stress_vec;
            results.Stress = von_mises; % Scalar for plotting
            results.Force = 0; 
        end

        % =================================================================
        % GetLocationMatrix (获取定位向量)
        % =================================================================
        function lm = GetLocationMatrix(obj)
            lm = zeros(12, 1);
            for i = 1:4
                % 每个节点只取前 3 个自由度 (u, v, w)
                % Node 1: lm(1:3) = BCode(1:3)
                % Node 2: lm(4:6) = BCode(1:3) ...
                
                idx_start = (i-1)*3 + 1;
                % 只提取 BCode 的 1-3 分量
                lm(idx_start : idx_start+2) = obj.Nodes(i).BCode(1:3);
            end
        end

        % =================================================================
        % CalcMass (计算质量矩阵)
        % =================================================================
        function m = CalcMass(obj)
            % 1. 获取坐标并计算体积
            x=zeros(4,1); y=zeros(4,1); z=zeros(4,1);
            for i=1:4, n=obj.Nodes(i); x(i)=n.XYZ(1); y(i)=n.XYZ(2); z(i)=n.XYZ(3); end
            V6 = det([1,x(1),y(1),z(1); 1,x(2),y(2),z(2); 1,x(3),y(3),z(3); 1,x(4),y(4),z(4)]);
            V = abs(V6)/6.0;

            % 2. 获取密度
            rho = 7850; 
            if isprop(obj.Material, 'Rho'), rho = obj.Material.Rho; end
            
            % 3. 一致质量矩阵
            factor = (rho * V) / 20.0;
            m = zeros(12, 12);
            for i = 1:4
                for j = 1:4
                    if i == j, val = 2; else, val = 1; end
                    % 填充 3x3 对角块
                    r = (i-1)*3+1; c = (j-1)*3+1;
                    m(r:r+2, c:c+2) = val * factor * eye(3);
                end
            end
        end
    end

    % =====================================================================
    % Private Helper Methods
    % =====================================================================
    methods (Access = private)
        function D = GetDMatrix(obj)
            E = obj.Material.E;
            NU = obj.Material.Nu; % Assuming property name is Nu or Poisson
            
            
            A1 = NU / (1 - NU);
            A2 = (1 - 2*NU) / 2 / (1 - NU);
            factor = E * (1 - NU) / (1 + NU) / (1 - 2*NU);
            
            D_raw = [ 1,  A1, A1,  0,  0,  0;
                      A1, 1,  A1,  0,  0,  0;
                      A1, A1, 1,   0,  0,  0;
                      0,  0,  0,   A2, 0,  0;
                      0,  0,  0,   0,  A2, 0;
                      0,  0,  0,   0,  0,  A2];
            D = D_raw * factor;
        end
    end
end