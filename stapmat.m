%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILE PATH: stapmat.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stapmat(inputFileName)
% STAPMAT 主程序 (支持静态与动力学分析)
%
% 用法:
%   stapmat('Data/stap90.in')

    if nargin < 1
        inputFileName = 'Data/verification_100.in';
    end

    % 设置路径 
    addpath(genpath('SRC')); 
    
    % 准备输出文件
    [pathstr, name, ~] = fileparts(inputFileName);
    outputFileName = fullfile(pathstr, [name, '.out']);
    pltFileName = fullfile(pathstr, [name, '.plt']);
    
    fidOut = fopen(outputFileName, 'w');
    if fidOut == -1, error('无法创建输出文件: %s', outputFileName); end
    finishup = onCleanup(@() fclose(fidOut));

    % --- 1. 读取数据 ---
    femDomain = Domain.Instance();
    if ~femDomain.ReadData(inputFileName)
        error('读取输入文件失败。');
    end

    % --- 2. 输出控制信息 ---
    fprintf(fidOut, ' C O N T R O L   I N F O R M A T I O N\n\n');
    fprintf(fidOut, '      NUMBER OF NODAL POINTS . . . . . . . (NUMNP)  = %10d\n', femDomain.NUMNP);
    fprintf(fidOut, '      NUMBER OF ELEMENT GROUPS . . . . . . (NUMEG)  = %10d\n', femDomain.NUMEG);
    fprintf(fidOut, '      NUMBER OF LOAD CASES . . . . . . . . (NLCASE) = %10d\n', femDomain.NLCASE);
    fprintf(fidOut, '      SOLUTION MODE  . . . . . . . . . . . (MODEX)  = %10d\n\n', femDomain.MODEX);

    % --- 3. 求解阶段 ---
    tic;
    
    % 3.1 组装刚度矩阵 (无论是静态还是动力学都需要)
    femDomain.AssembleStiffnessMatrix();
    
    % 根据输入文件选择求解器
    if femDomain.MODEX == 1
        
        if femDomain.AnalysisType == 1
            % --- 动力学分析 ---
            fprintf(' [Mode] Dynamic Analysis Selected (Input File Configured).\n');
            
            p = femDomain.DynParams; % 从 Domain 获取读取的参数
            
            % 实例化广义-alpha 动力学求解器
            solver = GeneralizedAlphaSolver(p.dt, p.nSteps, p.rho_inf, p.alpha, p.beta);

            % 将 Tecplot 文件名传给求解器
            solver.OutputFileName = pltFileName;
            
            % 执行求解
            solver.Solve(femDomain);
            
        else
            % --- 静力学分析 ---
            fprintf(' [Mode] Static Analysis Selected.\n');
            
            % 实例化静态求解器
            solver = StaticSolver();
            solver.Solve(femDomain);
            % 静力学计算完成后，立即输出 Tecplot 文件
            fprintf(' Writing Static Result to Tecplot: %s ...\n', pltFileName);
            WriteTecplotLine(femDomain, pltFileName, 0.0, true);% 参数: domain, 文件名, 时间0.0, 是否新文件(true)
        end
    end
    
        % 执行求解 (多态调用)
        solver.Solve(femDomain); 
    
    totalTime = toc;

    % --- 4. 输出最终结果 ---
    % 对于动力学，这里输出的是最后一帧的位移；对于静力学，是平衡位置
    fprintf(fidOut, '\n D I S P L A C E M E N T S (Final State)\n\n');
    fprintf(fidOut, '    NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT\n');
    for i = 1:femDomain.NUMNP
        node = femDomain.NodeList(i);
        d = node.Displacement;
        fprintf(fidOut, '%8d      %15.6e   %15.6e   %15.6e\n', ...
            node.ID, d(1), d(2), d(3));
    end

    fprintf(fidOut, '\n S O L U T I O N   T I M E   L O G (SEC)\n\n');
    fprintf(fidOut, '     TOTAL TIME . . . . . . . . . . . . . . . . . . = %10.2f\n', totalTime);
    
    fprintf('计算完成！结果已写入: %s\n', outputFileName);
end