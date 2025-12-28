% BEAMELEMENT 3D Beam Element Class
%
% Purpose:
%   Represent a 2-node 3D beam element with 6 DOFs per node.
%   Compute stiffness matrix considering bending, torsion, and axial deformation.
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
%
% Programmed by: 宋博
classdef BeamElement < Element
    properties
        ReferenceVector % Reference Vector (参考向量，用于定义局部 y-z 平面)
    end
    
    methods
        % Constructor (构造函数)
        function obj = BeamElement()
            obj@Element(); 
            obj.ReferenceVector = zeros(3, 1); % 初始化为零向量
        end
        
        % Read element data(读取单元数据)
        % 格式: ElementID Node1 Node2 MatID RefX RefY RefZ
        function Read(obj, fid, expectedID, materialList, globalNodeList)
            % 读取一行数据
            lineStr = fgetl(fid);
            data = str2num(lineStr); %#ok<ST2NM>
            
            % 1. Check ID (校验单元编号)
            inputID = round(data(1));
            if inputID ~= expectedID
                error('Element ID mismatch. Expected: %d, Read: %d', expectedID, inputID);
            end
            obj.ID = inputID;
            
            % 2. Set Nodes (设置节点连接关系)
            n1 = round(data(2));
            n2 = round(data(3));
            obj.SetNodes(globalNodeList, [n1, n2]);
            
            % 3. Set Material (设置材料属性)
            matID = round(data(4));
            obj.Material = materialList(matID); 
            
            % 4. Set Reference Vector (设置参考向量)
            if length(data) < 7
                error('Beam element %d data insufficient. Need Ref Vector (X,Y,Z).', inputID);
            end
            obj.ReferenceVector(1) = data(5);
            obj.ReferenceVector(2) = data(6);
            obj.ReferenceVector(3) = data(7);
        end
        
        % CalcStiffness (计算全局单元刚度矩阵)
        % Output: k (12x12 matrix)
        function k = CalcStiffness(obj)
            % 1. Get Local Stiffness & Rotation Matrix (获取局部矩阵和旋转矩阵)
            [k_loc, T, ~] = obj.GetLocalMatrices();
           
            % 2. Transform to Global (坐标转换: K = T' * k_loc * T)
            k = T' * k_loc * T;
        end
        
        % GetLocationMatrix (获取定位向量/方程号)
        % Output: lm (12x1 vector)
        function lm = GetLocationMatrix(obj)
            lm = zeros(12, 1);
            lm(1:6)  = obj.Nodes(1).BCode; % Node 1: DOFs 1-6
            lm(7:12) = obj.Nodes(2).BCode; % Node 2: DOFs 7-12
        end

        % CalcMass (计算一致质量矩阵)
        % Output: m (12x12 matrix)
        function m = CalcMass(obj)
            [~, T, L] = obj.GetLocalMatrices(); 
            % 读取材料信息
            mat = obj.Material;
            rho = mat.Density;
            A   = mat.Area;
            J   = mat.J;

            totalMass = rho * A * L;

            % Axial Terms (u)
            m_axial = (totalMass / 6.0) * [2, 1; 1, 2];     

            % Torsion Terms (theta_x) - Polar mass moment of inertia
            % J*rho*L / 6 * [2, 1; 1, 2]
            m_torsion = (rho * J * L / 6.0) * [2, 1; 1, 2];
            
            % Bending Coefficients (Bernoulli Beam)
            c = totalMass / 420.0;
            
            % Block for v-theta_z (XY Plane)
            m_bend = c * [156,     22*L,    54,     -13*L;
                          22*L,    4*L^2,   13*L,   -3*L^2;
                          54,      13*L,    156,    -22*L;
                          -13*L,   -3*L^2,  -22*L,  4*L^2];
            
            % 3. Assemble Local Mass Matrix (12x12)
            m_loc = zeros(12, 12);
            
            % Axial (1, 7)
            m_loc([1,7], [1,7]) = m_axial;
            
            % Torsion (4, 10)
            m_loc([4,10], [4,10]) = m_torsion;
            
            % Bending XY (v, th_z -> 2, 6, 8, 12)
            idx_xy = [2, 6, 8, 12];
            m_loc(idx_xy, idx_xy) = m_bend;
            
            % Bending XZ (w, th_y -> 3, 5, 9, 11)
            % For XZ plane, rotation signs are flipped for off-diagonals involving theta_y
            % standard consistent mass matrix for w-thy
            m_bend_xz = c * [156,     -22*L,   54,     13*L;
                             -22*L,   4*L^2,   -13*L,  -3*L^2;
                             54,      13*L,    156,    22*L;
                             13*L,    -3*L^2,  22*L,   4*L^2];
            idx_xz = [3, 5, 9, 11];
            m_loc(idx_xz, idx_xz) = m_bend_xz;
            
            % 4. Transform to Global
            m = T' * m_loc * T;
        end
               
        % CalcStress (计算单元内力和应力，输出最大正应力)
        % Output: results (结构体，包含 ForceVector 和 Stress)
         function results = CalcStress(obj, globalU)
            % 1. Get Matrices (获取矩阵)
            [k_loc, T, ~] = obj.GetLocalMatrices();
            
            % 2. 提取单元相关的 12 个全局位移
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
            
            % 3. 转换到局部坐标系: u_loc = T * u_g
            u_loc = T * u_g;
            
            % 4. 计算局部内力: f = k * u
            % f_loc 包含了两个节点的力: 
            % [Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2]'
            f_loc = k_loc * u_loc;
            
            % 5. 提取内力
            % 轴力 N (取节点2的力，拉为正)
            N  = f_loc(7); 
            
            % 弯矩 
            % 局部坐标系下：
            % 绕 Y 轴弯矩: M_y1 = f_loc(5),  M_y2 = f_loc(11)
            % 绕 Z 轴弯矩: M_z1 = f_loc(6),  M_z2 = f_loc(12)
            My_max = max(abs(f_loc(5)), abs(f_loc(11)));
            Mz_max = max(abs(f_loc(6)), abs(f_loc(12)));

            % 提取扭矩 (Torsion) - Mx
            % 理论上梁两端扭矩平衡但反向，取绝对值
            Torsion = abs(f_loc(10));

            % 提取合成弯矩 (Resultant Bending Moment)
            % 需要检查两端 (节点1和节点2)，取最大值
            % Node 1: M1 = sqrt(My1^2 + Mz1^2)
            M1 = sqrt(f_loc(5)^2 + f_loc(6)^2);
            % Node 2: M2 = sqrt(My2^2 + Mz2^2)
            M2 = sqrt(f_loc(11)^2 + f_loc(12)^2);
            
            MaxMoment = max(M1, M2);
            
            % 6. 计算最大正应力: Axial + Bending
            % Formula: sigma = |N/A| + |My * z_max / Iy| + |Mz * y_max / Iz|
            mat = obj.Material;
            A  = mat.Area;
            Iy = mat.Iy;
            Iz = mat.Iz;
            
            % --- 估算截面边缘距离 c ---
            % 基于"实心矩形"假设进行估算 (c = sqrt(3 * I/A))
            c_z = sqrt(3 * Iy / A); % 对应绕 Y 轴弯曲，Z 方向的边缘距离
            c_y = sqrt(3 * Iz / A); % 对应绕 Z 轴弯曲，Y 方向的边缘距离
            
            % 计算截面模量 W (Section Modulus)
            Wy = Iy / c_z;
            Wz = Iz / c_y;
            
            % 叠加应力 (绝对值叠加，代表截面上受力最大的点的应力)
            sigma_axial = abs(N / A);
            sigma_bend  = My_max / Wy + Mz_max / Wz;
            sigma_max = sigma_axial + sigma_bend;
            
            % 7. 打包结果
            results.ForceVector = [N, f_loc(8), f_loc(9), f_loc(10), f_loc(11), f_loc(12)];
            results.Force  = N;
            results.Stress = sigma_max; % 最大组合正应力
            results.Torsion = Torsion;
            results.Moment  = MaxMoment; % 输出扭矩和弯矩
        end
        
    end

    
    % =====================================================================
    % Private Helper Methods(得到局部单元刚度矩阵、转换矩阵、单元长度)
    % =====================================================================
    methods (Access = private)
        
        function [k_loc, T, L] = GetLocalMatrices(obj)
            % 1. Geometric Properties (几何属性)
            node1 = obj.Nodes(1);
            node2 = obj.Nodes(2);
            d = node2.XYZ - node1.XYZ;
            L = norm(d);
            
            if L <= 1e-12
                error('Element %d has zero length!', obj.ID);
            end
            
            % 2. Material Properties (材料属性)
            mat = obj.Material; 
            E = mat.E;   % Young's Modulus (杨氏模量)
            G = mat.G;   % Shear Modulus (剪切模量)
            A = mat.Area;% Area (截面积)
            Iy = mat.Iy; % Moment of Inertia Y (惯性矩 Iy)
            Iz = mat.Iz; % Moment of Inertia Z (惯性矩 Iz)
            J = mat.J;   % Torsional Constant (扭转常数)
            
            % 3. Stiffness Coefficients (刚度系数)
            X  = E * A / L;             % Axial stiffness
            Y1 = 12 * E * Iz / L^3;     % Bending Z terms
            Y2 = 6 * E * Iz / L^2;
            Y3 = 4 * E * Iz / L;
            Y4 = 2 * E * Iz / L;
            Z1 = 12 * E * Iy / L^3;     % Bending Y terms
            Z2 = 6 * E * Iy / L^2;
            Z3 = 4 * E * Iy / L;
            Z4 = 2 * E * Iy / L;
            S  = G * J / L;             % Torsion stiffness
            
            % 4. Fill Local Stiffness Matrix (填充局部刚度矩阵 12x12)
            k_loc = zeros(12, 12);
            
           % ... (前面计算 X, Y1, Y2, Z1, Z2 等系数的部分保持不变)

            % ---------------------------------------------------------
            % 修正后的刚度矩阵赋值 (Standard Right-Hand Rule)
            % ---------------------------------------------------------
            
            % [1] Axial (轴向)
            k_loc(1,1) = X;   k_loc(1,7) = -X;
            k_loc(7,1) = -X;  k_loc(7,7) = X;

            % [2] Torsion (扭转)
            k_loc(4,4) = S;   k_loc(4,10)= -S;
            k_loc(10,4)= -S;  k_loc(10,10)= S;

            % [3] Bending XY Plane (v - theta_z) -> 保持正号逻辑
            % Node 1 Internal
            k_loc(2,2) = Y1;
            k_loc(6,6) = Y3;
            k_loc(2,6) = Y2;  % Positive
            k_loc(6,2) = Y2;
            
            % Node 2 Internal
            k_loc(8,8) = Y1;
            k_loc(12,12)= Y3;
            k_loc(8,12) = -Y2; % Negative (Symmetric antisymmetry)
            k_loc(12,8) = -Y2;
            
            % Interaction
            k_loc(2,8) = -Y1; k_loc(8,2) = -Y1;
            k_loc(6,12)= Y4;  k_loc(12,6)= Y4;
            k_loc(2,12)= Y2;  k_loc(12,2)= Y2;
            k_loc(6,8) = -Y2; k_loc(8,6) = -Y2;

            % [4] Bending XZ Plane (w - theta_y) -> 关键修正：符号必须与 XY 平面相反
            % Node 1 Internal
            k_loc(3,3) = Z1;
            k_loc(5,5) = Z3;
            k_loc(3,5) = -Z2;  
            k_loc(5,3) = -Z2;  
            
            % Node 2 Internal
            k_loc(9,9) = Z1;
            k_loc(11,11)= Z3;
            k_loc(9,11) = Z2;  
            k_loc(11,9) = Z2;  
            
            % Interaction
            k_loc(3,9) = -Z1; k_loc(9,3) = -Z1;
            k_loc(5,11)= Z4;  k_loc(11,5)= Z4;
            
            k_loc(3,11)= -Z2; 
            k_loc(11,3)= -Z2; 
            k_loc(5,9) = Z2;  
            k_loc(9,5) = Z2;  
            
            % ---------------------------------------------------------
            
            % 5. Coordinate Transformation Matrix T (坐标转换矩阵)
            % Local x axis (u-axis)
            vec_x = d / L;
            
            % Reference vector (参考向量)
            vec_ref = obj.ReferenceVector;
            
            % Local z axis (w-axis) = x cross ref
            vec_z = cross(vec_x, vec_ref);
            
            % Handle Singularity: if ref is parallel to beam, use global Y
            if norm(vec_z) < 1e-10
                % Warning suppressed for standard workflow, but logic is handled
                if abs(vec_x(1)) < 0.9 && abs(vec_x(2)) > 0.9 % Beam is Y-axis
                    vec_z = cross(vec_x, [0; 0; 1]); % Use global Z
                else
                    vec_z = cross(vec_x, [0; 1; 0]); % Use global Y
                end
            end
            vec_z = vec_z / norm(vec_z);
            
            % Local y axis (v-axis) = z cross x
            vec_y = cross(vec_z, vec_x);
            vec_y = vec_y / norm(vec_y);
            
            % Build 3x3 Rotation Block
            T_sub = [vec_x'; vec_y'; vec_z'];
            
            % Assemble 12x12 Transformation Matrix T
            % T = blkdiag(T_sub, T_sub, T_sub, T_sub) equivalent
            T = zeros(12, 12);
            T(1:3, 1:3)     = T_sub;
            T(4:6, 4:6)     = T_sub;
            T(7:9, 7:9)     = T_sub;
            T(10:12, 10:12) = T_sub;

        end
    end
end