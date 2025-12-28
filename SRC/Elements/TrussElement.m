% 
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
% 修改：程云志

classdef TrussElement < Element    
    methods
        function obj = TrussElement()
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
            
            n1 = round(data(2));
            n2 = round(data(3));
            obj.SetNodes(globalNodeList, [n1, n2]);
            
            matID = round(data(4));
            obj.Material = materialList(matID); 
        end
        
        % 单元刚度矩阵，为了实现梁-杆混合，也改成12×12大小
        function k = CalcStiffness(obj)
            node1 = obj.Nodes(1);
            node2 = obj.Nodes(2);
            
            dx = node2.XYZ(1) - node1.XYZ(1);
            dy = node2.XYZ(2) - node1.XYZ(2);
            dz = node2.XYZ(3) - node1.XYZ(3);
            
            L2 = dx*dx + dy*dy + dz*dz;
            L = sqrt(L2); % 单元长度
            
            % 读取材料属性
            E = obj.Material.E;
            A = obj.Material.Area;
            
            % k = E*A/L
            coeff = E * A / L;
            
            % 转换矩阵T(3x1)
            nx = dx / L;
            ny = dy / L;
            nz = dz / L;
            T = [nx; ny; nz];
            
            % 局部 3x3 子块
            k_sub = coeff * (T * T'); 
            
            % 组装到 12x12 矩阵
            % 索引对应：Node1(1-3), Node2(7-9)。转动自由度(4-6, 10-12)为零。
            
            k = zeros(12, 12);
            
            % Node 1 - Node 1 (左上)
            k(1:3, 1:3) = k_sub;
            
            % Node 1 - Node 2 (右上)
            k(1:3, 7:9) = -k_sub;
            
            % Node 2 - Node 1 (左下)
            k(7:9, 1:3) = -k_sub;
            
            % Node 2 - Node 2 (右下)
            k(7:9, 7:9) = k_sub;
        end

        % =================================================================
        % CalcThermalLoad (计算桁架单元的热等效节点力)
        % =================================================================
        function f = CalcThermalLoad(obj)
            
            node1 = obj.Nodes(1);
            node2 = obj.Nodes(2);
            
            dx = node2.XYZ(1) - node1.XYZ(1);
            dy = node2.XYZ(2) - node1.XYZ(2);
            dz = node2.XYZ(3) - node1.XYZ(3);
            L = sqrt(dx^2 + dy^2 + dz^2);
            
            
            if L < 1e-12, L=1; end 
            n = [dx/L; dy/L; dz/L];
            
            
            T1 = 0; T2 = 0;
            if isprop(node1, 'Temperature'), T1 = node1.Temperature; end
            if isprop(node2, 'Temperature'), T2 = node2.Temperature; end
            
            T_avg = (T1 + T2) / 2.0;
            
            
            E = obj.Material.E;
            A = obj.Material.Area;
            Alpha = 0.0;
            if isprop(obj.Material, 'Alpha'), Alpha = obj.Material.Alpha; end
            
            F_mag = E * A * Alpha * T_avg;
            
            
            
            f1 = -F_mag * n;
            f2 =  F_mag * n;
            
            f = zeros(12, 1);
            f(1:3) = f1; 
            f(7:9) = f2; 
        end
        
        % 计算桁架单元的应力和内力
        function results = CalcStress(obj, globalU)
            node1 = obj.Nodes(1);
            node2 = obj.Nodes(2);
            
            dx = node2.XYZ(1) - node1.XYZ(1);
            dy = node2.XYZ(2) - node1.XYZ(2);
            dz = node2.XYZ(3) - node1.XYZ(3);
            L = sqrt(dx^2 + dy^2 + dz^2);
            
            % 2. 提取单元节点的全局位移 (只取平动 DOF 1-3)
            u1_g = [0;0;0]; u2_g = [0;0;0];
            
            for i=1:3
                if node1.BCode(i) > 0, u1_g(i) = globalU(node1.BCode(i)); end
                if node2.BCode(i) > 0, u2_g(i) = globalU(node2.BCode(i)); end
            end
            
            % 3. 计算局部坐标系下的轴向变形 
            T = [dx/L, dy/L, dz/L]; 
            delta_L = T * (u2_g - u1_g);

            strain_total = delta_L / L;
            
            % 4. 热应变 (Thermal Strain)
            Alpha = 0.0;
            if isprop(obj.Material, 'Alpha'), Alpha = obj.Material.Alpha; end
            
            T1 = 0; T2 = 0;
            if isprop(node1, 'Temperature'), T1 = node1.Temperature; end
            if isprop(node2, 'Temperature'), T2 = node2.Temperature; end
            T_avg = (T1 + T2) / 2.0;
            
            strain_thermal = Alpha * T_avg;
            
            % 5. 计算机械应力: Sigma = E * (Eps_total - Eps_thermal)
            E = obj.Material.E;
            stress = E * (strain_total - strain_thermal);
            
            force = stress * obj.Material.Area;
            
            % 6. 打包结果
            results.Force = force;   
            results.Stress = stress; 
            % 新增兼容接口：扭矩和弯矩设为0
            results.Torsion = 0.0;
            results.Moment  = 0.0;
        end
        
        function lm = GetLocationMatrix(obj)
            
            lm = zeros(12, 1);
            lm(1:6)  = obj.Nodes(1).BCode;
            lm(7:12) = obj.Nodes(2).BCode;
        end

        % 质量矩阵计算 （集中质量）
        function m = CalcMass(obj)
            node1 = obj.Nodes(1); 
            node2 = obj.Nodes(2);
            L = norm(node2.XYZ - node1.XYZ);
            rho = obj.Material.Density; 
            A = obj.Material.Area;
            
            % 总质量
            totalMass = rho * A * L;
            
            % 集中质量(仅分配给平动自由度DOF1,2,3&7,8,9)
            m_lumped = totalMass / 2.0;
            m = zeros(12, 12);
            % Node 1 (u,v,w)
            m(1,1) = m_lumped; m(2,2) = m_lumped; m(3,3) = m_lumped;
            % Node 2 (u,v,w)
            m(7,7) = m_lumped; m(8,8) = m_lumped; m(9,9) = m_lumped;
        end
        
    end
end