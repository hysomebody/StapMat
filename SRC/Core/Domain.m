% DOMAIN Root class for the Finite Element Problem
%
% 功能：
%   核心管理类。
%   1. 管理全局数据：节点列表 (NodeList)、单元组 (ElemGroups)、材料集 (MaterialSets)、载荷工况 (LoadCases)。
%   2. 调度分析流程：读取数据 -> 编号 -> 组装矩阵 -> 求解 -> 回写结果。 
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
        Title       % 标题行
        MODEX       % 0: Data Check, 1: Execution
        NLCASE      % Number of Load Cases 载荷工况数量
        AnalysisType % 0: Static, 1: Dynamic
        DynParams    % Struct: .dt, .nSteps, .rho_inf, .alpha, .beta 动力学参数结构体
        
        % Infrastructure
        NodeList    % 节点对象数组(存储所有 Node)
        ElemGroups  % 单元组，每个 Cell 存储一种类型的单元数组
        MaterialSets% 材料集
        
        % Solution Data
        NEQ         % 总方程数 (系统自由度总和)
        GlobalK     % 全局刚度矩阵 (稀疏矩阵格式)
        GlobalM     % 全局质量矩阵 (稀疏矩阵格式)
        LoadCases   % 载荷工况数组
        
        % Counters
        NUMNP       % Number of Nodal Points 节点总数
        NUMEG       % Number of Element Groups 单元组总数
    end
    
    methods (Static)
        function obj = Instance()
            persistent uniqueInstance
            if isempty(uniqueInstance)
                uniqueInstance = Domain();
            end
            obj = uniqueInstance;
        end
    end
    
    methods (Access = private)
        function obj = Domain()
            obj.NodeList = Node.empty;
            obj.ElemGroups = {};
            obj.MaterialSets = {};
        end
    end
    
    methods (Access = public)

        % 每次读取数据前强制清除所有旧数据
        function Reset(obj)
            obj.NodeList = [];
            obj.LoadCases = [];
            obj.GlobalK = [];
            obj.GlobalM = [];
            obj.NEQ = 0;
            obj.NLCASE = 0;
            obj.NUMNP = 0;
            obj.NUMEG = 0;
            obj.ElemGroups = [] ;
        end

        function success = ReadData(obj, filename)
            obj.Reset();
            fid = fopen(filename, 'r');
            if fid == -1
                error('Cannot open file: %s', filename);
            end
            
            try         
                fprintf('Input phase ...\n');
                obj.Title = fgetl(fid);
                lineStr = fgetl(fid);
                tmp = str2num(lineStr); %#ok<ST2NM>
                obj.NUMNP  = round(tmp(1));
                obj.NUMEG  = round(tmp(2));
                obj.NLCASE = round(tmp(3));
                obj.MODEX  = round(tmp(4));

                % 这行参数大于5说明是动力学问题求解
                if length(tmp) >= 5
                    obj.AnalysisType = round(tmp(5));
                else
                    obj.AnalysisType = 0; % 默认为静力学
                end
                
                % 如果是动力学，读取下一行参数
                if obj.AnalysisType == 1
                    lineStr = fgetl(fid);
                    dtmp = str2num(lineStr);
                    if isempty(dtmp)
                        error('Dynamic Analysis selected but parameters line is missing.');
                    end
                    
                    obj.DynParams.dt      = dtmp(1);
                    obj.DynParams.nSteps  = round(dtmp(2));
                    obj.DynParams.rho_inf = dtmp(3);
                    
                    % 读取阻尼参数
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

                % 读取节点坐标与边界条件
                obj.ReadNodalPoints(fid);
                
                % 计算全局方程号
                obj.CalculateEquationNumber();
                
                % 读取载荷工况
                obj.ReadLoadCases(fid);
                
                % 读取单元组 (材料与单元连接)
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
        
        function ReadNodalPoints(obj, fid)
            fprintf('Reading %d Nodal Points...\n', obj.NUMNP);
            obj.NodeList = Node.empty(0, obj.NUMNP); % Pre-allocate
            
            for i = 1:obj.NUMNP
                obj.NodeList(i) = Node();
                obj.NodeList(i).Read(fid, i);
            end
        end
        
        % 计算方程号
        function CalculateEquationNumber(obj)
            fprintf('Calculating Equation Numbers (NDF=6)...\n');
            obj.NEQ = 0;
            
            for i = 1:obj.NUMNP
                node = obj.NodeList(i);
                for dof = 1:Node.NDF % 6个自由度
                    if node.BCode(dof) == 0
                        obj.NEQ = obj.NEQ + 1; % BC = 0 ,自由，分配新的方程号
                        node.BCode(dof) = obj.NEQ;
                    else
                        node.BCode(dof) = 0;
                    end
                end
            end
            fprintf('Total Active Equations (NEQ): %d\n', obj.NEQ);
        end
        
        % 读取载荷工况
        function ReadLoadCases(obj, fid)
            
            obj.LoadCases = struct('Nodes', [], 'DOFs', [], 'Mags', []);
            
            for i = 1:obj.NLCASE
                tmp = [];
                while isempty(tmp)
                    lineStr = fgetl(fid);
                    if ~ischar(lineStr) % 遇到文件尾
                        error('Unexpected End of File while reading Load Case %d header.', i);
                    end
                    tmp = str2num(lineStr); 
                end

                lcNum = tmp(1);
                numLoads = tmp(2);
                
                fprintf('Reading Load Case %d with %d loads...\n', lcNum, numLoads);
                
                nodes = zeros(numLoads, 1);
                dofs = zeros(numLoads, 1);
                mags = zeros(numLoads, 1);
                
                for j = 1:numLoads
                    lData = [];
                    while isempty(lData)
                        lLine = fgetl(fid);
                        if ~ischar(lLine)
                            error('Unexpected End of File while reading Load data (Case %d, Load %d).', i, j);
                        end
                        lData = str2num(lLine); 
                    end
                    % ------------------------------------------------
                    nodes(j) = lData(1);
                    dofs(j) = lData(2);
                    mags(j) = lData(3);
                end
                
                obj.LoadCases(i).Nodes = nodes;
                obj.LoadCases(i).DOFs = dofs;
                obj.LoadCases(i).Mags = mags;
            end
        end
        
        % 读单元组
        function ReadElements(obj, fid)
            fprintf('Reading %d Element Groups...\n', obj.NUMEG);
            
            for grp = 1:obj.NUMEG
                tmp = [];
                while isempty(tmp)
                    lineStr = fgetl(fid);
                    if ~ischar(lineStr)
                         error('Unexpected End of File reading Element Group %d.', grp);
                    end
                    tmp = str2num(lineStr); 
                end
                
                eType = tmp(1);
                nEle  = tmp(2);
                nMat  = tmp(3);
                
                fprintf('Group %d: Type=%d, Elements=%d, Materials=%d\n', grp, eType, nEle, nMat);
                
                if eType == 1 % Truss Element
                    mats = TrussMaterial.empty(0, nMat);
                    for m = 1:nMat
                        mats(m) = TrussMaterial();
                        mats(m).Read(fid, m);
                    end
                    obj.MaterialSets{grp} = mats;
                    
                    eles = TrussElement.empty(0, nEle);
                    for e = 1:nEle
                        eles(e) = TrussElement();
                        eles(e).Read(fid, e, mats, obj.NodeList);
                    end
                    obj.ElemGroups{grp} = eles;
                    
                elseif eType == 2 %  Beam Element
                    mats = BeamMaterial.empty(0, nMat);
                    for m = 1:nMat
                        mats(m) = BeamMaterial();
                        mats(m).Read(fid, m);
                    end
                    obj.MaterialSets{grp} = mats;
                    
                    eles = BeamElement.empty(0, nEle);
                    for e = 1:nEle
                        eles(e) = BeamElement();
                        eles(e).Read(fid, e, mats, obj.NodeList);
                    end
                    obj.ElemGroups{grp} = eles;

                elseif eType == 3 % Tetra Element 
                    mats = TetraMaterial.empty(0, nMat);
                    for m = 1:nMat
                        mats(m) = TetraMaterial();
                        mats(m).Read(fid, m);
                    end
                    obj.MaterialSets{grp} = mats;
                    
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
        
        % 组装全局刚度矩阵
        function AssembleStiffnessMatrix(obj)
            fprintf('Assembling Global Stiffness Matrix...\n');
            
            totalEle = 0;
            for g = 1:length(obj.ElemGroups)
                totalEle = totalEle + length(obj.ElemGroups{g});
            end
            estimatedNZ = totalEle * 144; 
            
            I_idx = zeros(estimatedNZ, 1);
            J_idx = zeros(estimatedNZ, 1);
            K_val = zeros(estimatedNZ, 1);
            
            count = 0;
            

            for g = 1:length(obj.ElemGroups)
                elements = obj.ElemGroups{g};
                
                for e = 1:length(elements)
                    elem = elements(e);
                    % 单元刚度矩阵
                    ke = elem.CalcStiffness();
                    % 定位向量
                    lm = elem.GetLocationMatrix();
                    for r = 1:12
                        rowEq = lm(r);
                        if rowEq == 0; continue; end 
                        
                        for c = 1:12
                            colEq = lm(c);
                            if colEq == 0; continue; end
                            
                            count = count + 1;
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
            
            if count == 0
                warning('Stiffness Matrix is empty!');
                obj.GlobalK = sparse(obj.NEQ, obj.NEQ);
                return;
            end
            I_idx = I_idx(1:count);
            J_idx = J_idx(1:count);
            K_val = K_val(1:count);
            
            max_I = max(I_idx);
            max_J = max(J_idx);
            max_Idx = max(max_I, max_J);
            
            if max_Idx > obj.NEQ
                fprintf('Error: Max Index (%d) > NEQ (%d). Resizing Global Matrix.\n', max_Idx, obj.NEQ);
                obj.NEQ = max_Idx; 
            elseif obj.NEQ > max_Idx
            end
            
            obj.GlobalK = sparse(I_idx, J_idx, K_val, obj.NEQ, obj.NEQ);
            
            fprintf('Assembly Done. Matrix Size: %dx%d, Non-zeros: %d\n', ...
                obj.NEQ, obj.NEQ, nnz(obj.GlobalK));
        end
        
        % 组装全局质量矩阵（动力学计算用）
        function AssembleMassMatrix(obj)
            fprintf('Assembling Global Mass Matrix...\n');
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
                    
                    % 单元质量矩阵
                    me = elem.CalcMass(); 
                    % 定位向量
                    lm = elem.GetLocationMatrix();

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
            
            I_idx = I_idx(1:count);
            J_idx = J_idx(1:count);
            M_val = M_val(1:count);
            
            obj.GlobalM = sparse(I_idx, J_idx, M_val, obj.NEQ, obj.NEQ);
            
            fprintf('Mass Assembly Done. Non-zeros: %d\n', nnz(obj.GlobalM));
        end
       
        % 组装全局力向量 F
        function [F_total, F_mech, F_thermal] = AssembleForce(obj, loadCaseIdx)
            if loadCaseIdx > obj.NLCASE
                error('Load Case Index out of range');
            end
            
            fprintf('Assembling Force Vector for Load Case %d (Separating Mech & Thermal)...\n', loadCaseIdx);
            
            F_mech = zeros(obj.NEQ, 1);
            F_thermal = zeros(obj.NEQ, 1);
            
            % --- 1. 施加集中载荷
            lc = obj.LoadCases(loadCaseIdx);
            nLoads = length(lc.Nodes);
            for i = 1:nLoads
                nodeID = lc.Nodes(i);
                dofDir = lc.DOFs(i); 
                mag    = lc.Mags(i);
                
                node = obj.NodeList(nodeID);
                eqNum = node.BCode(dofDir);

                if eqNum > 0
                    F_mech(eqNum) = F_mech(eqNum) + mag;
                end
            end
            
            % --- 2. 施加 热载荷
            for g = 1:length(obj.ElemGroups)
                elements = obj.ElemGroups{g};
                if isempty(elements), continue; end
                
                firstElem = elements(1);
                if ismethod(firstElem, 'CalcThermalLoad')
                    
                    for e = 1:length(elements)
                        elem = elements(e);
                        
                        f_ele = elem.CalcThermalLoad();
                        
                        lm = elem.GetLocationMatrix();
                        
                        for i = 1:length(lm)
                            eq = lm(i);
                            if eq > 0
                                F_thermal(eq) = F_thermal(eq) + f_ele(i);
                            end
                        end
                    end
                end
            end
            
            F_total = F_mech + F_thermal;
            
            fprintf('   Force Assembly Done. Mech Norm: %.2e, Thermal Norm: %.2e\n', norm(F_mech), norm(F_thermal));
        end
        
        % 将计算出的全局位移向量 U_val 分发回各个节点
        function UpdateNodalDisplacements(obj, U_val)
            % fprintf('Updating Nodal Displacements...\n');
            % fprintf('   NODE          X-DISP          Y-DISP          Z-DISP\n');
            for i = 1:obj.NUMNP
                node = obj.NodeList(i);
                dispVec = zeros(6, 1); % 6 DOF
                
                for dof = 1:6
                    eq = node.BCode(dof);
                    if eq > 0
                        dispVec(dof) = U_val(eq);
                    end
                end
                node.Displacement = dispVec;   
                
            % % 打印节点平动分量，验证结果
            % fprintf(' %6d  %14.6e  %14.6e  %14.6e\n', ...
            %     node.ID, dispVec(1), dispVec(2), dispVec(3));

            end
        end

    end
end