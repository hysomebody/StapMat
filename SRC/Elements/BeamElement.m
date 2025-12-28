% 
% 功能：3D 欧拉-伯努利梁（Euler-Bernoulli Beam）单元。
%   每个节点包含 6 个自由度 (X, Y, Z 平动 + RX, RY, RZ 转动)。
%   支持双向弯曲、扭转和轴向变形耦合。
%
% 核心逻辑：
%   1. 根据截面属性 (Area, Iy, Iz, J, E, G) 计算局部刚度矩阵 k_loc。
%   2. 利用参考向量 (Reference Vector) 构建坐标转换矩阵 T。
%   3. 转换得到全局刚度矩阵 K = T' * k_loc * T。
%   4. 根据全局位移反算单元内力 (CalcStress)。
%
% Call procedures:
%   ./Element.m - Element() (Base Constructor)
%   ./Element.m - SetNodes()
%   ./Node.m
%   ./BeamMaterial.m
%
% Called by:
%   ./Domain.m

classdef BeamElement < Element
    properties
        ReferenceVector % 参考向量，用于定义局部 y-z 平面
    end
    
    methods
        function obj = BeamElement()
            obj@Element(); 
            obj.ReferenceVector = zeros(3, 1); % 初始化为零向量
        end
        
        function Read(obj, fid, expectedID, materialList, globalNodeList)
            % 读取一行数据
            lineStr = fgetl(fid);
            data = str2num(lineStr); %#ok<ST2NM>
            
            inputID = round(data(1));
            if inputID ~= expectedID
                error('Element ID mismatch. Expected: %d, Read: %d', expectedID, inputID);
            end
            obj.ID = inputID;
            
            n1 = round(data(2));
            n2 = round(data(3));
            obj.SetNodes(globalNodeList, [n1, n2]);
            
            matID = round(data(4));
            obj.Material = materialList(matID); 
            
            % 设置参考向量
            if length(data) < 7
                error('Beam element %d data insufficient. Need Ref Vector (X,Y,Z).', inputID);
            end
            obj.ReferenceVector(1) = data(5);
            obj.ReferenceVector(2) = data(6);
            obj.ReferenceVector(3) = data(7);
        end
        
        % CalcStiffness (计算全局单元刚度矩阵)
        function k = CalcStiffness(obj)
            [k_loc, T, ~] = obj.GetLocalMatrices();
           
            % 坐标转换: K = T' * k_loc * T
            k = T' * k_loc * T;
        end
        
        % 获取定位向量/方程号
        % lm (12x1)
        function lm = GetLocationMatrix(obj)
            lm = zeros(12, 1);
            lm(1:6)  = obj.Nodes(1).BCode; % Node 1: DOFs 1-6
            lm(7:12) = obj.Nodes(2).BCode; % Node 2: DOFs 7-12
        end

        % CalcMass (计算一致质量矩阵)
        %  m (12x12)
        function m = CalcMass(obj)
            [~, T, L] = obj.GetLocalMatrices(); 
            mat = obj.Material;
            rho = mat.Density;
            A   = mat.Area;
            J   = mat.J;

            totalMass = rho * A * L;

            % 轴向分量
            m_axial = (totalMass / 6.0) * [2, 1; 1, 2];     

            % 扭转分量
            m_torsion = (rho * J * L / 6.0) * [2, 1; 1, 2];
            
            % 弯曲分量
            c = totalMass / 420.0;
            
            % XY 平面弯曲 (对应 v 位移和 theta_z 转角)
            m_bend = c * [156,     22*L,    54,     -13*L;
                          22*L,    4*L^2,   13*L,   -3*L^2;
                          54,      13*L,    156,    -22*L;
                          -13*L,   -3*L^2,  -22*L,  4*L^2];

            % XZ 平面弯曲(对应 w 位移和 theta_y 转角)
            m_bend_xz = c * [156,     -22*L,   54,     13*L;
                             -22*L,   4*L^2,   -13*L,  -3*L^2;
                             54,      13*L,    156,    22*L;
                             13*L,    -3*L^2,  22*L,   4*L^2];
            
            % 局部质量矩阵
            m_loc = zeros(12, 12);

            m_loc([1,7], [1,7]) = m_axial;

            m_loc([4,10], [4,10]) = m_torsion;

            idx_xy = [2, 6, 8, 12];
            m_loc(idx_xy, idx_xy) = m_bend;

            idx_xz = [3, 5, 9, 11];
            m_loc(idx_xz, idx_xz) = m_bend_xz;
            
            m = T' * m_loc * T;
        end
               
        % CalcStress (计算单元内力和应力，输出最大正应力)
         function results = CalcStress(obj, globalU)
            [k_loc, T, ~] = obj.GetLocalMatrices();
            
            u_g = zeros(12, 1);
            % Node 1 DOFs (1-6)
            for i = 1:6
                id = obj.Nodes(1).BCode(i);
                if id > 0, u_g(i) = globalU(id); end
            end
            % Node 2 DOFs (7-12)
            for i = 1:6
                id = obj.Nodes(2).BCode(i);
                if id > 0, u_g(i+6) = globalU(id); end
            end
            
            % 转换到局部坐标系: u_loc = T * u_g
            u_loc = T * u_g;
            
            % 计算局部内力: f = k * u
            % f_loc 包含了两个节点的力: 
            % [Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2]'
            f_loc = k_loc * u_loc;
            
            % 提取内力
            % 轴力 N (取节点2的力，拉为正)
            N  = f_loc(7); 
            
            % 弯矩 (取两端绝对值的最大值)
            My_max = max(abs(f_loc(5)), abs(f_loc(11)));
            Mz_max = max(abs(f_loc(6)), abs(f_loc(12)));

            % 扭矩
            Torsion = abs(f_loc(10));

            % 合成弯矩
            % Node 1: M1 = sqrt(My1^2 + Mz1^2)
            M1 = sqrt(f_loc(5)^2 + f_loc(6)^2);
            % Node 2: M2 = sqrt(My2^2 + Mz2^2)
            M2 = sqrt(f_loc(11)^2 + f_loc(12)^2);
            
            MaxMoment = max(M1, M2);
            
            % 计算最大组合正应力
            % sigma = |N/A| + |My * z_max / Iy| + |Mz * y_max / Iz|
            mat = obj.Material;
            A  = mat.Area;
            Iy = mat.Iy;
            Iz = mat.Iz;
            
            c_z = sqrt(3 * Iy / A); 
            c_y = sqrt(3 * Iz / A); 
            
            % 截面模量 W
            Wy = Iy / c_z;
            Wz = Iz / c_y;
            
            % 叠加应力 (绝对值叠加，代表截面上受力最大的点的应力)
            sigma_axial = N / A;
            sigma_bend  = My_max / Wy + Mz_max / Wz;
            sigma_1 = sigma_axial + sigma_bend;
            sigma_2 = sigma_axial - sigma_bend;
            if abs(sigma_1) > abs(sigma_2)
                sigma_max = sigma_1;
            else
                sigma_max = sigma_2;
            end
            % sigma_max = sigma_axial + sigma_bend;
            
            results.ForceVector = [N, f_loc(8), f_loc(9), f_loc(10), f_loc(11), f_loc(12)];
            results.Force  = N;
            results.Stress = sigma_max; 
            results.Torsion = Torsion;
            results.Moment  = MaxMoment; 
        end
        
    end


    methods (Access = private)
        
        function [k_loc, T, L] = GetLocalMatrices(obj)
            node1 = obj.Nodes(1);
            node2 = obj.Nodes(2);
            d = node2.XYZ - node1.XYZ;
            L = norm(d);
            
            if L <= 1e-12
                error('Element %d has zero length!', obj.ID);
            end
            

            mat = obj.Material; 
            E = mat.E;   
            G = mat.G;   
            A = mat.Area;
            Iy = mat.Iy; 
            Iz = mat.Iz; 
            J = mat.J;   
            
            
            X  = E * A / L;             
            Y1 = 12 * E * Iz / L^3;     
            Y2 = 6 * E * Iz / L^2;
            Y3 = 4 * E * Iz / L;
            Y4 = 2 * E * Iz / L;
            Z1 = 12 * E * Iy / L^3;     
            Z2 = 6 * E * Iy / L^2;
            Z3 = 4 * E * Iy / L;
            Z4 = 2 * E * Iy / L;
            S  = G * J / L;             
            
            k_loc = zeros(12, 12);
            
            k_loc(1,1) = X;   k_loc(1,7) = -X;
            k_loc(7,1) = -X;  k_loc(7,7) = X;

            k_loc(4,4) = S;   k_loc(4,10)= -S;
            k_loc(10,4)= -S;  k_loc(10,10)= S;

            k_loc(2,2) = Y1;
            k_loc(6,6) = Y3;
            k_loc(2,6) = Y2;  
            k_loc(6,2) = Y2;
            
            k_loc(8,8) = Y1;
            k_loc(12,12)= Y3;
            k_loc(8,12) = -Y2; 
            k_loc(12,8) = -Y2;
            
            k_loc(2,8) = -Y1; k_loc(8,2) = -Y1;
            k_loc(6,12)= Y4;  k_loc(12,6)= Y4;
            k_loc(2,12)= Y2;  k_loc(12,2)= Y2;
            k_loc(6,8) = -Y2; k_loc(8,6) = -Y2;

            k_loc(3,3) = Z1;
            k_loc(5,5) = Z3;
            k_loc(3,5) = -Z2;  
            k_loc(5,3) = -Z2;  
            
            k_loc(9,9) = Z1;
            k_loc(11,11)= Z3;
            k_loc(9,11) = Z2;  
            k_loc(11,9) = Z2;  
            
            k_loc(3,9) = -Z1; k_loc(9,3) = -Z1;
            k_loc(5,11)= Z4;  k_loc(11,5)= Z4;
            
            k_loc(3,11)= -Z2; 
            k_loc(11,3)= -Z2; 
            k_loc(5,9) = Z2;  
            k_loc(9,5) = Z2;  
            
            vec_x = d / L;
            
            vec_ref = obj.ReferenceVector;
            
            vec_z = cross(vec_x, vec_ref);
            
            if norm(vec_z) < 1e-10
                if abs(vec_x(1)) < 0.9 && abs(vec_x(2)) > 0.9 
                    vec_z = cross(vec_x, [0; 0; 1]); 
                else
                    vec_z = cross(vec_x, [0; 1; 0]); 
                end
            end
            vec_z = vec_z / norm(vec_z);
            
            vec_y = cross(vec_z, vec_x);
            vec_y = vec_y / norm(vec_y);

            T_sub = [vec_x'; vec_y'; vec_z'];
            
            T = zeros(12, 12);
            T(1:3, 1:3)     = T_sub;
            T(4:6, 4:6)     = T_sub;
            T(7:9, 7:9)     = T_sub;
            T(10:12, 10:12) = T_sub;

        end
    end
end