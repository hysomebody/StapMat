% TRUSSELEMENT Truss (Bar) Element Class
%
% Purpose:
%   Represent a 2-node truss element.
%   Compute stiffness matrix for 3D truss.
% 功能：实现桁架单元的刚度矩阵计算
% 
% Call procedures:
%   ./Element.m - Element()
%   ./Element.m - SetNodes()
%   ./Node.m
%   ./TrussMaterial.m
%
% Called by:
%   ./Domain.m

classdef TrussElement < Element    
    methods
        % Constructor
        function obj = TrussElement()
            obj@Element(); % Call base constructor
        end
        
        % Read element data
        % Format: ElementID Node1 Node2 MaterialID
        function Read(obj, fid, expectedID, materialList, globalNodeList)
            % called by: Domain (parsing element block)
            
            lineStr = fgetl(fid);
            data = str2num(lineStr); %#ok<ST2NM>
            
            % 1. Check ID
            inputID = round(data(1));
            if inputID ~= expectedID
                error('Element ID mismatch. Expected: %d, Read: %d', expectedID, inputID);
            end
            obj.ID = inputID;
            
            % 2. Set Nodes
            n1 = round(data(2));
            n2 = round(data(3));
            obj.SetNodes(globalNodeList, [n1, n2]);
            
            % 3. Set Material
            matID = round(data(4));
            % Assuming materialList is an array of Material objects
            % and matID matches the index (1-based)
            obj.Material = materialList(matID); 
        end
        
        % Calculate 12x12 Stiffness Matrix (NDOF=6 compatible)
        function k = CalcStiffness(obj)
            % called by: Domain.AssembleStiffnessMatrix
            
            % 1. Get Geometry
            node1 = obj.Nodes(1);
            node2 = obj.Nodes(2);
            
            dx = node2.XYZ(1) - node1.XYZ(1);
            dy = node2.XYZ(2) - node1.XYZ(2);
            dz = node2.XYZ(3) - node1.XYZ(3);
            
            L2 = dx*dx + dy*dy + dz*dz;
            L = sqrt(L2);
            
            % 2. Get Material Properties
            % Cast to TrussMaterial to access Area (though not strictly needed in dynamic languages like MATLAB, it's good practice)
            E = obj.Material.E;
            A = obj.Material.Area;
            
            % 3. Compute Basic Stiffness Coefficient
            % k = E*A/L
            coeff = E * A / L;
            
            % 4. Direction Cosines / Transformation vector
            % Original code logic (TrussStiff.m):
            % ST(1) = dx/L^2 ... This effectively computes (dx/L)*(1/L)
            % Let's stick to standard FEM formulation: n = [dx, dy, dz] / L
            nx = dx / L;
            ny = dy / L;
            nz = dz / L;
            
            % T vector for one node (3x1)
            T = [nx; ny; nz];
            
            % 5. Build Local 3x3 Block (k_sub = coeff * T * T')
            k_sub = coeff * (T * T'); 
            
            % 6. Assemble into 12x12 Matrix (for 2 nodes * 6 DOF)
            % Indices for Translation DOFs: 
            % Node 1: 1, 2, 3 (X, Y, Z)
            % Node 2: 7, 8, 9 (X, Y, Z)
            % Rotational DOFs (4-6, 10-12) are zero for Truss.
            
            k = zeros(12, 12);
            
            % Node 1 - Node 1 (Top-Left)
            k(1:3, 1:3) = k_sub;
            
            % Node 1 - Node 2 (Top-Right)
            k(1:3, 7:9) = -k_sub;
            
            % Node 2 - Node 1 (Bottom-Left)
            k(7:9, 1:3) = -k_sub;
            
            % Node 2 - Node 2 (Bottom-Right)
            k(7:9, 7:9) = k_sub;
        end
        
        % 计算桁架单元的应力和内力
        function results = CalcStress(obj, globalU)
            % 1. 获取节点和几何信息
            node1 = obj.Nodes(1);
            node2 = obj.Nodes(2);
            
            dx = node2.XYZ(1) - node1.XYZ(1);
            dy = node2.XYZ(2) - node1.XYZ(2);
            dz = node2.XYZ(3) - node1.XYZ(3);
            L = sqrt(dx^2 + dy^2 + dz^2);
            
            % 2. 提取单元节点的全局位移 (只取平动 DOF 1-3)
            % Node 1 (Equation numbers check is handled by extraction from full vector)
            u1_g = [0;0;0]; u2_g = [0;0;0];
            
            for i=1:3
                if node1.BCode(i) > 0, u1_g(i) = globalU(node1.BCode(i)); end
                if node2.BCode(i) > 0, u2_g(i) = globalU(node2.BCode(i)); end
            end
            
            % 3. 计算局部坐标系下的轴向变形 (Transform to Local)
            % Transformation vector T = [nx, ny, nz]
            T = [dx/L, dy/L, dz/L]; 
            
            % local deformation = T * (u2 - u1)
            delta_L = T * (u2_g - u1_g);
            
            % 4. 计算应变和应力
            strain = delta_L / L;
            E = obj.Material.E;
            stress = E * strain;
            force = stress * obj.Material.Area;
            
            % 5. 打包结果
            results.Force = force;   
            results.Stress = stress; 
            % 新增兼容接口：扭矩和弯矩设为0
            results.Torsion = 0.0;
            results.Moment  = 0.0;
        end
        
        % Return the Location Matrix (Equation numbers)
        function lm = GetLocationMatrix(obj)
            % called by: Domain.AssembleStiffnessMatrix
            
            % Gather BCode from all nodes
            % Node 1 DOFs: 1-6 -> mapped to lm(1:6)
            % Node 2 DOFs: 1-6 -> mapped to lm(7:12)
            
            lm = zeros(12, 1);
            lm(1:6)  = obj.Nodes(1).BCode;
            lm(7:12) = obj.Nodes(2).BCode;
        end

        % 实现质量矩阵计算
        function m = CalcMass(obj)
            % 1. 获取单元长度、密度、截面积等材料信息
            node1 = obj.Nodes(1); 
            node2 = obj.Nodes(2);
            L = norm(node2.XYZ - node1.XYZ);
            
            rho = obj.Material.Density; 
            A = obj.Material.Area;
            
            % 2. 总质量
            totalMass = rho * A * L;
            
            % 3. 集中质量(仅分配给平动自由度DOF1,2,3&7,8,9)
            m_lumped = totalMass / 2.0;
            
            m = zeros(12, 12);
            % Node 1 (u,v,w)
            m(1,1) = m_lumped; m(2,2) = m_lumped; m(3,3) = m_lumped;
            % Node 2 (u,v,w)
            m(7,7) = m_lumped; m(8,8) = m_lumped; m(9,9) = m_lumped;
        end
        
    end
end