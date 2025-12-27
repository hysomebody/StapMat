% DOMAIN Root class for the Finite Element Problem
%
% Purpose:
%   Manage global data (NodeList, Elements, Materials, LoadCases).
%   Orchestrate the analysis process (Read -> Number -> Assemble -> Solve).
% 功能：整个程序的大管家，读取材料、读取单元、计算方程号和组装刚度矩阵
% 
% Call procedures:
%   ./Node.m - Read()
%   ./TrussMaterial.m - Read()
%   ./TrussElement.m - Read()
%   ./TrussElement.m - CalcStiffness()
%   ./TrussElement.m - GetLocationMatrix()
%
% Called by:
%   stapmat.m (Main Script) - To initialize and read data

classdef Domain < handle

    properties
        % Control Data
        Title       % Problem Title
        MODEX       % 0: Data Check, 1: Execution
        NLCASE      % Number of Load Cases 载荷数量
        AnalysisType % 0: Static, 1: Dynamic
        DynParams    % Struct: .dt, .nSteps, .rho_inf, .alpha, .beta
        
        % Infrastructure
        NodeList    % (Array of Node) All nodes
        ElemGroups  % (Cell array) Stores arrays of Elements (e.g., {TrussEleArray, BeamEleArray})
        MaterialSets% (Cell array) Stores arrays of Materials
        
        % Solution Data
        NEQ         % Total Number of Equations
        GlobalK     % Global Stiffness Matrix (Sparse)
        GlobalM     % Global Mass Matrix (Sparse)
        LoadCases   % (Struct Array) Placeholder for load data
        
        % Counters
        NUMNP       % Number of Nodal Points
        NUMEG       % Number of Element Groups
    end
    
    methods (Static)
        % Singleton Instance Access
        function obj = Instance()
            persistent uniqueInstance
            if isempty(uniqueInstance)
                uniqueInstance = Domain();
            end
            obj = uniqueInstance;
        end
    end
    
    methods (Access = private)
        % Private Constructor
        function obj = Domain()
            obj.NodeList = Node.empty;
            obj.ElemGroups = {};
            obj.MaterialSets = {};
        end
    end
    
    methods (Access = public)
        
        % Main entry point to read input file
        function success = ReadData(obj, filename)
            fid = fopen(filename, 'r');
            if fid == -1
                error('Cannot open file: %s', filename);
            end
            
            try
                % 取新文件前，清除上一次分析的残留矩阵
                obj.GlobalK = [];
                obj.GlobalM = [];
          
                fprintf('Input phase ...\n');
                
                % 1. Read Heading
                obj.Title = fgetl(fid);
                
                % 2. Read Control Line
                lineStr = fgetl(fid);
                tmp = str2num(lineStr); %#ok<ST2NM>
                obj.NUMNP  = round(tmp(1));
                obj.NUMEG  = round(tmp(2));
                obj.NLCASE = round(tmp(3));
                obj.MODEX  = round(tmp(4));

                if length(tmp) >= 5
                    obj.AnalysisType = round(tmp(5));
                else
                    obj.AnalysisType = 0; % 默认为静力学
                end
                
                % 如果是动力学，读取下一行参数
                if obj.AnalysisType == 1
                    lineStr = fgetl(fid);
                    dtmp = str2num(lineStr); %#ok<ST2NM>
                    if isempty(dtmp)
                        error('Dynamic Analysis selected but parameters line is missing.');
                    end
                    
                    obj.DynParams.dt      = dtmp(1);
                    obj.DynParams.nSteps  = round(dtmp(2));
                    obj.DynParams.rho_inf = dtmp(3);
                    
                    % 读取阻尼 (可选)
                    if length(dtmp) >= 5
                        obj.DynParams.alpha = dtmp(4);
                        obj.DynParams.beta  = dtmp(5);
                    else
                        obj.DynParams.alpha = 0.0;
                        obj.DynParams.beta  = 0.0;
                    end
                    
                    fprintf('Dynamic Params Read: dt=%.2e, N=%d, rho=%.2f\n', ...
                        obj.DynParams.dt, obj.DynParams.nSteps, obj.DynParams.rho_inf);
                end

                % 3. Read Nodal Points
                obj.ReadNodalPoints(fid);
                
                % 4. Calculate Equation Numbers (Crucial for 6DOF)
                obj.CalculateEquationNumber();
                
                % 5. Read Load Cases
                obj.ReadLoadCases(fid);
                
                % 6. Read Element Groups
                obj.ReadElements(fid);
                
                fclose(fid);
                success = true;
                
            catch ME
                if fid ~= -1
                    fclose(fid);
                end
                rethrow(ME);
            end
        end
        
        % Read Node data block
        function ReadNodalPoints(obj, fid)
            fprintf('Reading %d Nodal Points...\n', obj.NUMNP);
            obj.NodeList = Node.empty(0, obj.NUMNP); % Pre-allocate
            
            for i = 1:obj.NUMNP
                obj.NodeList(i) = Node();
                obj.NodeList(i).Read(fid, i);
            end
        end
        
        % Calculate global equation numbers based on BCs
        function CalculateEquationNumber(obj)
            fprintf('Calculating Equation Numbers (NDF=6)...\n');
            obj.NEQ = 0;
            
            % Loop over all nodes
            for i = 1:obj.NUMNP
                node = obj.NodeList(i);
                for dof = 1:Node.NDF % Loop 1 to 6
                    if node.BCode(dof) == 0
                        % If BC is 0 (Free), assign a new equation number
                        obj.NEQ = obj.NEQ + 1;
                        node.BCode(dof) = obj.NEQ;
                    else
                        % If BC is 1 (Fixed), assign 0
                        node.BCode(dof) = 0;
                    end
                end
            end
            fprintf('Total Active Equations (NEQ): %d\n', obj.NEQ);
        end
        
        % Read Load Cases (Simplified implementation for now)
        function ReadLoadCases(obj, fid)
            % Format based on stap90.in:
            % Line 1: LoadCaseNum, NumConcentratedLoads
            % Following lines: NodeID, DOF, Magnitude
            
            obj.LoadCases = struct('Nodes', [], 'DOFs', [], 'Mags', []);
            
            for i = 1:obj.NLCASE
                lineStr = fgetl(fid);
                tmp = str2num(lineStr); %#ok<ST2NM>
                lcNum = tmp(1);
                numLoads = tmp(2);
                
                fprintf('Reading Load Case %d with %d loads...\n', lcNum, numLoads);
                
                nodes = zeros(numLoads, 1);
                dofs = zeros(numLoads, 1);
                mags = zeros(numLoads, 1);
                
                for j = 1:numLoads
                    lLine = fgetl(fid);
                    lData = str2num(lLine); %#ok<ST2NM>
                    nodes(j) = lData(1);
                    dofs(j) = lData(2);
                    mags(j) = lData(3);
                end
                
                obj.LoadCases(i).Nodes = nodes;
                obj.LoadCases(i).DOFs = dofs;
                obj.LoadCases(i).Mags = mags;
            end
        end
        
        % Read Element Groups
        function ReadElements(obj, fid)
            fprintf('Reading %d Element Groups...\n', obj.NUMEG);
            
            % Loop over each group
            for grp = 1:obj.NUMEG
                % Read Group Control Line: Type NumEle NumMat
                lineStr = fgetl(fid);
                tmp = str2num(lineStr); %#ok<ST2NM>
                if isempty(tmp); lineStr = fgetl(fid); tmp = str2num(lineStr); end % Skip empty lines if any
                
                eType = tmp(1);
                nEle  = tmp(2);
                nMat  = tmp(3);
                
                fprintf('Group %d: Type=%d, Elements=%d, Materials=%d\n', grp, eType, nEle, nMat);
                
                if eType == 1 % Truss Element
                    % 1. Read Materials for this group
                    mats = TrussMaterial.empty(0, nMat);
                    for m = 1:nMat
                        mats(m) = TrussMaterial();
                        mats(m).Read(fid, m);
                    end
                    obj.MaterialSets{grp} = mats;
                    
                    % 2. Read Elements for this group
                    eles = TrussElement.empty(0, nEle);
                    for e = 1:nEle
                        eles(e) = TrussElement();
                        % Pass global NodeList and local MaterialList
                        eles(e).Read(fid, e, mats, obj.NodeList);
                    end
                    obj.ElemGroups{grp} = eles;
                    
                elseif eType == 2 %  Beam Element
                    % 1. Read Materials for this group (BeamMaterial)
                    mats = BeamMaterial.empty(0, nMat);
                    for m = 1:nMat
                        mats(m) = BeamMaterial();
                        mats(m).Read(fid, m);
                    end
                    obj.MaterialSets{grp} = mats;
                    
                    % 2. Read Elements for this group (BeamElement)
                    eles = BeamElement.empty(0, nEle);
                    for e = 1:nEle
                        eles(e) = BeamElement();
                        % Pass global NodeList (to link nodes) and local MaterialList
                        eles(e).Read(fid, e, mats, obj.NodeList);
                    end
                    obj.ElemGroups{grp} = eles;

                elseif eType == 3 % Tetra Element (New Implementation)
                    % 1. Read Materials (TetraMaterial)
                    mats = TetraMaterial.empty(0, nMat);
                    for m = 1:nMat
                        mats(m) = TetraMaterial();
                        mats(m).Read(fid, m);
                    end
                    obj.MaterialSets{grp} = mats;
                    
                    % 2. Read Elements (TetraElement)
                    eles = TetraElement.empty(0, nEle);
                    for e = 1:nEle
                        eles(e) = TetraElement();
                        eles(e).Read(fid, e, mats, obj.NodeList);
                    end
                    obj.ElemGroups{grp} = eles;

                else
                    error('Unknown Element Type: %d. (Supported: 1=Truss, 2=Beam, 3=Tetra)', eType);
                end
            end
        end
        
        % 组装全局刚度矩阵 Assemble Global Stiffness Matrix
        function AssembleStiffnessMatrix(obj)
            fprintf('Assembling Global Stiffness Matrix...\n');
            
            % Estimate non-zeros to pre-allocate sparse triplet
            % Conservative guess: 12x12 per element
            totalEle = 0;
            for g = 1:length(obj.ElemGroups)
                totalEle = totalEle + length(obj.ElemGroups{g});
            end
            estimatedNZ = totalEle * 144; 
            
            I_idx = zeros(estimatedNZ, 1);
            J_idx = zeros(estimatedNZ, 1);
            K_val = zeros(estimatedNZ, 1);
            
            count = 0;
            
            % Loop over all groups
            for g = 1:length(obj.ElemGroups)
                elements = obj.ElemGroups{g};
                
                % Loop over elements in group
                for e = 1:length(elements)
                    elem = elements(e);
                    
                    % 1. Get Element Stiffness (12x12)
                    ke = elem.CalcStiffness();
                    
                    % 2. Get Location Matrix (Equation numbers)
                    lm = elem.GetLocationMatrix();
                    
                    % 3. Assemble into triplet vectors
                    for r = 1:12
                        rowEq = lm(r);
                        if rowEq == 0; continue; end % Constrained DOF
                        
                        for c = 1:12
                            colEq = lm(c);
                            if colEq == 0; continue; end
                            
                            count = count + 1;
                            % 动态扩容防止溢出 (Safety check)
                            if count > length(I_idx)
                                I_idx(end*2) = 0;
                                J_idx(end*2) = 0;
                                K_val(end*2) = 0;
                            end

                            I_idx(count) = rowEq;
                            J_idx(count) = colEq;
                            K_val(count) = ke(r, c);
                        end
                    end
                end
            end
            
            % Trim unused pre-allocation
            if count == 0
                warning('Stiffness Matrix is empty!');
                obj.GlobalK = sparse(obj.NEQ, obj.NEQ);
                return;
            end
            I_idx = I_idx(1:count);
            J_idx = J_idx(1:count);
            K_val = K_val(1:count);
            
            % --- 诊断与自动修正 ---
            max_I = max(I_idx);
            max_J = max(J_idx);
            max_Idx = max(max_I, max_J);
            
            if max_Idx > obj.NEQ
                fprintf('Error: Max Index (%d) > NEQ (%d). Resizing Global Matrix.\n', max_Idx, obj.NEQ);
                % 强制修正维度，避免报错，允许程序继续运行以查看结果
                obj.NEQ = max_Idx; 
            elseif obj.NEQ > max_Idx
                % 这通常是正常的（有些方程可能没有任何刚度贡献，尽管少见），但值得注意
                % fprintf('Info: NEQ (%d) > Max Index in K (%d).\n', obj.NEQ, max_Idx);
            end
            % --------------------
            
            % Create Sparse Matrix
            obj.GlobalK = sparse(I_idx, J_idx, K_val, obj.NEQ, obj.NEQ);
            
            fprintf('Assembly Done. Matrix Size: %dx%d, Non-zeros: %d\n', ...
                obj.NEQ, obj.NEQ, nnz(obj.GlobalK));
        end
        
        % 组装全局质量矩阵 (类似 AssembleStiffnessMatrix)
        function AssembleMassMatrix(obj)
            fprintf('Assembling Global Mass Matrix...\n');
            % 预估非零元 (跟刚度矩阵一样)
            totalEle = 0;
            for g = 1:length(obj.ElemGroups)
                totalEle = totalEle + length(obj.ElemGroups{g});
            end
            estimatedNZ = totalEle * 144; 
            
            I_idx = zeros(estimatedNZ, 1);
            J_idx = zeros(estimatedNZ, 1);
            M_val = zeros(estimatedNZ, 1);
            
            count = 0;
            
            for g = 1:length(obj.ElemGroups)
                elements = obj.ElemGroups{g};
                for e = 1:length(elements)
                    elem = elements(e);
                    
                    % 1. 获取单元质量矩阵
                    me = elem.CalcMass(); % 调用单元的接口
                    
                    % 2. 获取定位向量
                    lm = elem.GetLocationMatrix();
                    
                    % 3. 组装
                    for r = 1:12
                        rowEq = lm(r);
                        if rowEq == 0; continue; end 
                        
                        for c = 1:12
                            colEq = lm(c);
                            if colEq == 0; continue; end
                            
                            count = count + 1;
                            I_idx(count) = rowEq;
                            J_idx(count) = colEq;
                            M_val(count) = me(r, c);
                        end
                    end
                end
            end
            
            % 裁剪并生成稀疏矩阵
            I_idx = I_idx(1:count);
            J_idx = J_idx(1:count);
            M_val = M_val(1:count);
            
            obj.GlobalM = sparse(I_idx, J_idx, M_val, obj.NEQ, obj.NEQ);
            
            fprintf('Mass Assembly Done. Non-zeros: %d\n', nnz(obj.GlobalM));
        end
       
        % 组装全局力向量 F
        function F = AssembleForce(obj, loadCaseIdx)
            if loadCaseIdx > obj.NLCASE
                error('Load Case Index out of range');
            end
            
            fprintf('Assembling Force Vector for Load Case %d...\n', loadCaseIdx);
            
            % 初始化全局力向量
            F = zeros(obj.NEQ, 1);
            
            % 获取该工况的载荷数据
            lc = obj.LoadCases(loadCaseIdx);
            nLoads = length(lc.Nodes);
            
            for i = 1:nLoads
                nodeID = lc.Nodes(i);
                dofDir = lc.DOFs(i); % 1=X, 2=Y, 3=Z...
                mag    = lc.Mags(i);
                
                % 找到该节点
                node = obj.NodeList(nodeID);
                
                % 获取该自由度对应的方程号
                eqNum = node.BCode(dofDir);
                
                % 如果方程号 > 0，说明是自由自由度，需要施加力
                if eqNum > 0
                    F(eqNum) = F(eqNum) + mag;
                else
                    fprintf('Warning: Load applied to constrained DOF (Node %d, DOF %d)\n', nodeID, dofDir);
                end
            end

            % --- 2. 施加 体积力/热载荷 (Thermal Loads) ---
            % 遍历所有单元组
            for g = 1:length(obj.ElemGroups)
                elements = obj.ElemGroups{g};
                if isempty(elements), continue; end
                
                % 检查该组单元是否支持热载荷计算 (即是否有 CalcThermalLoad 方法)
                % 这样可以兼容 Truss/Beam 单元（如果它们还没写这个方法的话）
                firstElem = elements(1);
                if ismethod(firstElem, 'CalcThermalLoad')
                    
                    for e = 1:length(elements)
                        elem = elements(e);
                        
                        % A. 计算单元热载荷向量 (12x1)
                        f_ele = elem.CalcThermalLoad();
                        
                        % B. 获取定位向量 (Location Matrix)
                        lm = elem.GetLocationMatrix();
                        
                        % C. 组装到全局向量 F
                        for i = 1:length(lm)
                            eq = lm(i);
                            if eq > 0
                                F(eq) = F(eq) + f_ele(i);
                            end
                        end
                    end
                end
            end
        end  
        
        % 将计算出的全局位移向量 U_val 分发回各个节点
        function UpdateNodalDisplacements(obj, U_val)
            fprintf('Updating Nodal Displacements...\n');
            fprintf('   NODE          X-DISP          Y-DISP          Z-DISP\n');
            
            for i = 1:obj.NUMNP
                node = obj.NodeList(i);
                dispVec = zeros(6, 1); % 6 DOF
                
                for dof = 1:6
                    eq = node.BCode(dof);
                    if eq > 0
                        dispVec(dof) = U_val(eq);
                    end
                end
            % 将计算结果存入节点对象 
                node.Displacement = dispVec;   
                
           %  % 打印节点平动分量，验证结果
           %  fprintf(' %6d  %14.6e  %14.6e  %14.6e\n', ...
           %      node.ID, dispVec(1), dispVec(2), dispVec(3));
            end
        end

    end
end