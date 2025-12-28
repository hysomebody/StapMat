% TETRAELEMENT 4-Node Tetrahedral Element Class (Linear)
%
% Purpose:
%   Represent a 4-node 3D solid element (Tetrahedron).
%   Compute stiffness matrix for 3D elasticity.
% 功能：实现四节点四面体单元（常应变单元 CST）。
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

        function k = CalcStiffness(obj)
            x = zeros(4,1); y = zeros(4,1); z = zeros(4,1);
            for i = 1:4
                node = obj.Nodes(i);
                x(i) = node.XYZ(1);
                y(i) = node.XYZ(2);
                z(i) = node.XYZ(3);
            end

            M_vol = [1, x(1), y(1), z(1);
                     1, x(2), y(2), z(2);
                     1, x(3), y(3), z(3);
                     1, x(4), y(4), z(4)];
            V6 = det(M_vol);
            V = abs(V6) / 6.0;

            B = zeros(6, 12);
            
            for i = 1:4
                j = mod(i, 4) + 1;
                m = mod(i + 1, 4) + 1;
                p = mod(i + 2, 4) + 1;

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

            D = obj.GetDMatrix();

            k = B' * D * B * V;
            
        end

        function results = CalcStress(obj, globalU)
            u_elem = zeros(12, 1);
            for i = 1:4
                node = obj.Nodes(i);
                for dof = 1:3
                    eqID = node.BCode(dof);
                    if eqID > 0
                        u_elem((i-1)*3 + dof) = globalU(eqID);
                    end
                end
            end

            x = zeros(4,1); y = zeros(4,1); z = zeros(4,1);
            for i = 1:4
                n = obj.Nodes(i); x(i)=n.XYZ(1); y(i)=n.XYZ(2); z(i)=n.XYZ(3);
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

            % Calculate Stress Vector
            D = obj.GetDMatrix();
            % Stress = D * Strain = D * B * u
            stress_vec = D * B * u_elem; 
            % stress_vec format: [sig_x, sig_y, sig_z, tau_yz, tau_zx, tau_xy]'

            % Calculate Von Mises Stress (Original Logic from Te4Stress)
            % I1 = sig1 + sig2 + sig3
            % I2 ... I3 ...
            % Using invariant formula directly from components:
            sx = stress_vec(1); sy = stress_vec(2); sz = stress_vec(3);
            t_yz = stress_vec(4); t_zx = stress_vec(5); t_xy = stress_vec(6);

            % Von Mises Formula:
            % sqrt(0.5 * [(sx-sy)^2 + (sy-sz)^2 + (sz-sx)^2 + 6*(txy^2 + tyz^2 + tzx^2)])
            vm_sq = 0.5 * ((sx-sy)^2 + (sy-sz)^2 + (sz-sx)^2 + 6*(t_xy^2 + t_yz^2 + t_zx^2));
            von_mises = sqrt(vm_sq);

            results.StressVector = stress_vec;
            results.Stress = von_mises; % Scalar for plotting
            results.Force = 0; % Force concept is vague for solid, maybe return 0
        end

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

        % CalcMass (计算质量矩阵 - 集中质量)
        function m = CalcMass(obj)
            % Volume
            x=zeros(4,1); y=zeros(4,1); z=zeros(4,1);
            for i=1:4, n=obj.Nodes(i); x(i)=n.XYZ(1); y(i)=n.XYZ(2); z(i)=n.XYZ(3); end
            V6 = det([1,x(1),y(1),z(1); 1,x(2),y(2),z(2); 1,x(3),y(3),z(3); 1,x(4),y(4),z(4)]);
            V = abs(V6)/6.0;

            rho = 7850; 
            if isprop(obj.Material, 'Rho'), rho = obj.Material.Rho; end
            totalMass = rho * V;
            m_node = totalMass / 4.0;
            
            % 12x12 质量矩阵 (仅平动)
            m = eye(12) * m_node;
        end
    end

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