% stapmat.m
function stapmat(inputFileName)
% STAPMAT 主程序
% 功能：
%   1. 初始化环境，设置路径。
%   2. 调用 Domain 读取有限元模型数据 (.in 文件)。
%   3. 格式化输出模型输入数据 (Echo Input)。
%   4. 组装全局刚度矩阵，求解位移。
%   5. 计算单元应力/内力。
%   6. 格式化输出求解结果 (.out 文件)。
%
% 用法:
%   stapmat('Data/stap90.in')
    
    % --- 1. 环境初始化 ---
    if nargin < 1
        inputFileName = 'Data/stap90.in'; % 默认输入文件
    end

    % 设置类文件搜索路径 (SRC 及其子文件夹)
    addpath(genpath('SRC'));
    
    % 准备输出文件
    [pathstr, name, ~] = fileparts(inputFileName);
    outputFileName = fullfile(pathstr, [name, '.out']);
    
    fidOut = fopen(outputFileName, 'w');
    if fidOut == -1, error('无法创建输出文件: %s', outputFileName); end
    % 使用 onCleanup 确保函数退出时自动关闭文件
    finishup = onCleanup(@() fclose(fidOut));
    
    % --- 2. 读取数据 (Input Phase) ---
    t_start = clock;
    
    femDomain = Domain.Instance();
    % 清理旧数据（如果是单例模式，需重置，此处假设每次运行重新读取）
    if ~femDomain.ReadData(inputFileName)
        error('读取输入文件失败。');
    end
    
    t_input = etime(clock, t_start);

    % --- 3. 输出模型信息 (Echo Input Data) ---
    
    % 3.1 控制信息 (Control Information)
    fprintf(fidOut, '\n C O N T R O L   I N F O R M A T I O N\n\n');
    fprintf(fidOut, '      NUMBER OF NODAL POINTS . . . . . . . . . . (NUMNP)  = %10d\n', femDomain.NUMNP);
    fprintf(fidOut, '      NUMBER OF ELEMENT GROUPS . . . . . . . . . (NUMEG)  = %10d\n', femDomain.NUMEG);
    fprintf(fidOut, '      NUMBER OF LOAD CASES . . . . . . . . . . . (NLCASE) = %10d\n', femDomain.NLCASE);
    fprintf(fidOut, '      SOLUTION MODE  . . . . . . . . . . . . . . (MODEX)  = %10d\n', femDomain.MODEX);
    fprintf(fidOut, '         EQ.0, DATA CHECK\n');
    fprintf(fidOut, '         EQ.1, EXECUTION\n\n');

    % 3.2 节点数据 (Nodal Point Data)
    fprintf(fidOut, ' N O D A L   P O I N T   D A T A\n\n');
    fprintf(fidOut, '       NODE          BOUNDARY                         NODAL POINT\n');
    fprintf(fidOut, '      NUMBER     CONDITION  CODES                     COORDINATES\n');
    fprintf(fidOut, '                    X    Y    Z               X            Y            Z\n');
    
    for i = 1:femDomain.NUMNP
        node = femDomain.NodeList(i);
        % 将方程号转换为 BC 代码显示 (0=Free, 1=Fixed)
        % 注意：Domain 中 BCode=0 是 Fixed, >0 是方程号(Free)。
        % 但原始 STAP90 输入/输出约定是：0=Free, 1=Fixed。
        % 我们需要根据 BCode 反推显示的 codes。
        
        displayCodes = zeros(1,3);
        for k=1:3
            if node.BCode(k) > 0 
                displayCodes(k) = 0; % 有方程号 -> Free
            else
                displayCodes(k) = 1; % 无方程号 -> Fixed
            end
        end
        
        fprintf(fidOut, '%10d      %5d%5d%5d      %13.3f%13.3f%13.3f\n', ...
            node.ID, displayCodes(1), displayCodes(2), displayCodes(3), ...
            node.XYZ(1), node.XYZ(2), node.XYZ(3));
    end

    % 3.3 方程号 (Equation Numbers)
    fprintf(fidOut, '\n\n EQUATION NUMBERS\n\n');
    fprintf(fidOut, '        NODE         DEGREES OF FREEDOM\n');
    fprintf(fidOut, '       NUMBER\n');
    fprintf(fidOut, '          N                  X         Y         Z\n');
    
    for i = 1:femDomain.NUMNP
        node = femDomain.NodeList(i);
        fprintf(fidOut, ' %10d         %10d%10d%10d\n', ...
            node.ID, node.BCode(1), node.BCode(2), node.BCode(3));
    end

    % 3.4 载荷工况数据 (Load Case Data)
    % 仅打印第一个工况作为示例（原程序逻辑）或遍历打印
    fprintf(fidOut, '\n\n L O A D   C A S E   D A T A\n');
    
    for lc = 1:femDomain.NLCASE
        loadCase = femDomain.LoadCases(lc);
        nLoads = length(loadCase.Nodes);
        
        fprintf(fidOut, '\n     LOAD CASE NUMBER . . . . . . . = %10d\n', lc);
        fprintf(fidOut, '     NUMBER OF CONCENTRATED LOADS . = %10d\n\n', nLoads);
        fprintf(fidOut, '        NODE       DIRECTION      LOAD\n');
        fprintf(fidOut, '       NUMBER                   MAGNITUDE\n');
        
        for k = 1:nLoads
            fprintf(fidOut, '%10d         %4d       %12.5e\n', ...
                loadCase.Nodes(k), loadCase.DOFs(k), loadCase.Mags(k));
        end
    end

    % 3.5 单元组数据 (Element Group Data)
    fprintf(fidOut, '\n\n E L E M E N T   G R O U P   D A T A\n');
    
    for grp = 1:femDomain.NUMEG
        elements = femDomain.ElemGroups{grp};
        materials = femDomain.MaterialSets{grp};
        
        if isempty(elements), continue; end
        firstElem = elements(1);
        
        % 判断单元类型代码
        elemType = 0;
        if isa(firstElem, 'TrussElement'), elemType = 1;
        elseif isa(firstElem, 'BeamElement'), elemType = 3; % STAP++ Beam is 3
        end
        
        fprintf(fidOut, '\n\n E L E M E N T   D E F I N I T I O N\n');
        fprintf(fidOut, '\n ELEMENT TYPE  . . . . . . . . . . . . .( NPAR(1) ) . . = %10d\n', elemType);
        fprintf(fidOut, '     EQ.1, TRUSS ELEMENTS\n');
        fprintf(fidOut, '     EQ.2, PLANE STRESS ELEMENTS\n');
        fprintf(fidOut, '     EQ.3, BEAM ELEMENTS\n');
        fprintf(fidOut, ' NUMBER OF ELEMENTS. . . . . . . . . . .( NPAR(2) ) . . = %10d\n', length(elements));
        
        % 打印材料属性
        fprintf(fidOut, '\n M A T E R I A L   D E F I N I T I O N\n');
        fprintf(fidOut, '\n NUMBER OF DIFFERENT SETS OF MATERIAL\n');
        fprintf(fidOut, ' AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . = %10d\n', length(materials));
        
        if elemType == 1 % Truss
            fprintf(fidOut, '  SET       YOUNG''S     CROSS-SECTIONAL\n');
            fprintf(fidOut, ' NUMBER     MODULUS          AREA\n');
            fprintf(fidOut, '               E              A\n');
            for m = 1:length(materials)
                mat = materials(m);
                fprintf(fidOut, '%5d    %12.5e  %14.6e\n', mat.ID, mat.E, mat.Area);
            end
        elseif elemType == 3 % Beam
            fprintf(fidOut, '  SET       YOUNG''S       SHEAR        CROSS-SECT      INERTIA-Y      INERTIA-Z\n');
            fprintf(fidOut, ' NUMBER     MODULUS       MODULUS         AREA            Iy             Iz\n');
            for m = 1:length(materials)
                mat = materials(m);
                fprintf(fidOut, '%5d    %10.3e   %10.3e   %10.3e   %10.3e   %10.3e\n', ...
                    mat.ID, mat.E, mat.G, mat.Area, mat.Iy, mat.Iz);
            end
        end
        
        % 打印单元连接信息
        fprintf(fidOut, '\n\n E L E M E N T   I N F O R M A T I O N\n');
        fprintf(fidOut, '\n      ELEMENT          NODE          NODE       MATERIAL\n');
        fprintf(fidOut, '      NUMBER-N           I             J       SET NUMBER\n');
        
        for e = 1:length(elements)
            el = elements(e);
            fprintf(fidOut, '%10d      %10d    %10d       %5d\n', ...
                el.ID, el.Nodes(1).ID, el.Nodes(2).ID, el.Material.ID);
        end
    end

    % --- 4. 求解阶段 (Solution Phase) ---
    t_before_assembly = clock;
    
    % 4.1 组装刚度矩阵
    femDomain.AssembleStiffnessMatrix();
    t_assembly = etime(clock, t_before_assembly);
    
    % 计算系统统计数据 (基于稀疏矩阵)
    K = femDomain.GlobalK;
    NEQ = femDomain.NEQ;
    NWK = nnz(K); % Number of non-zero elements
    
    % 计算带宽
    [i_idx, j_idx] = find(K);
    bandwidths = abs(i_idx - j_idx);
    MK = max(bandwidths); % Maximum half bandwidth
    % Mean bandwidth 估算
    MM = round(NWK / NEQ); 

    fprintf(fidOut, '\n\n  TOTAL SYSTEM DATA\n\n');
    fprintf(fidOut, '     NUMBER OF EQUATIONS . . . . . . . . . . . . . .(NEQ) = %10d\n', NEQ);
    fprintf(fidOut, '     NUMBER OF MATRIX ELEMENTS . . . . . . . . . . .(NWK) = %10d\n', NWK);
    fprintf(fidOut, '     MAXIMUM HALF BANDWIDTH  . . . . . . . . . . . .(MK ) = %10d\n', MK);
    fprintf(fidOut, '     MEAN HALF BANDWIDTH . . . . . . . . . . . . . .(MM ) = %10d\n', MM);

    % 4.2 求解位移 (Solve)
    t_before_solve = clock;
    
    % 为简单起见，这里演示 Load Case 1 的求解
    % 如果有多个工况，应该是一个循环
    if femDomain.MODEX == 1
        solver = StaticSolver();
        
        % 假设只处理第一个工况
        currentLC = 1;
        F = femDomain.AssembleForce(currentLC);
        
        % 处理稀疏矩阵
        U = K \ F;
        
        % 更新节点位移
        femDomain.UpdateNodalDisplacements(U);
    end
    
    t_solve = etime(clock, t_before_solve);

    % --- 5. 输出结果 (Result Output) ---
    
    % 5.1 输出位移
    fprintf(fidOut, '\n\n LOAD CASE %3d\n', 1);
    fprintf(fidOut, '\n D I S P L A C E M E N T S\n\n');
    fprintf(fidOut, '       NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT\n');
    
    for i = 1:femDomain.NUMNP
        node = femDomain.NodeList(i);
        d = node.Displacement;
        % 仅打印前三个平动自由度
        fprintf(fidOut, ' %10d        %18.6e%18.6e%18.6e\n', ...
            node.ID, d(1), d(2), d(3));
    end

    % 5.2 输出单元应力/内力
    t_before_stress = clock;
    
    for grp = 1:femDomain.NUMEG
        elements = femDomain.ElemGroups{grp};
        if isempty(elements), continue; end
        
        fprintf(fidOut, '\n\n  S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P %4d\n\n', grp);
        
        % 根据单元类型调整表头
        if isa(elements(1), 'TrussElement')
            fprintf(fidOut, '       ELEMENT             FORCE            STRESS\n');
            fprintf(fidOut, '       NUMBER\n');
            
            for e = 1:length(elements)
                el = elements(e);
                res = el.CalcStress(U); % 计算内力
                fprintf(fidOut, ' %10d           %13.6e    %13.6e\n', ...
                    el.ID, res.Force, res.Stress);
            end
            
        elseif isa(elements(1), 'BeamElement')
            fprintf(fidOut, '       ELEMENT             AXIAL            SHEAR-Y        SHEAR-Z        TORSION        MOMENT-Y       MOMENT-Z\n');
            fprintf(fidOut, '       NUMBER              FORCE\n');
            
            for e = 1:length(elements)
                el = elements(e);
                res = el.CalcStress(U);
                f = res.ForceVector; % [N, Qy, Qz, Mx, My, Mz]
                fprintf(fidOut, ' %10d        %12.4e   %12.4e   %12.4e   %12.4e   %12.4e   %12.4e\n', ...
                    el.ID, f(1), f(2), f(3), f(4), f(5), f(6));
            end
        end
    end
    
    t_stress = etime(clock, t_before_stress);
    total_time = etime(clock, t_start);

    % --- 6. 时间日志 (Time Log) ---
    fprintf(fidOut, '\n\n S O L U T I O N   T I M E   L O G   I N   S E C\n\n');
    fprintf(fidOut, '     TIME FOR INPUT PHASE  . . . . . . . . . . . . . . = %12.2f\n', t_input);
    fprintf(fidOut, '     TIME FOR CALCULATION OF STIFFNESS MATRIX  . . . . = %12.2f\n', t_assembly);
    fprintf(fidOut, '     TIME FOR FACTORIZATION AND LOAD CASE SOLUTIONS  . = %12.2f\n', t_solve);
    fprintf(fidOut, '     TIME FOR STRESS CALCULATIONS  . . . . . . . . . . = %12.2f\n\n', t_stress);
    fprintf(fidOut, '      T O T A L   S O L U T I O N   T I M E  . . . . . = %12.2f\n', total_time);

    WriteTecplot(femDomain)
    
    % 在控制台打印计算完成的消息
    fprintf('计算完成！结果已写入: %s\n', outputFileName);
end