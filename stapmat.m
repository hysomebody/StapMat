function stapmat(inputFileName)
% STAPMAT 主程序 (支持静态与动力学分析)
%
% 用法:
%   命令行键入stapmat('Data/stap90.in')，或其他输入文件名（确保输入文件要放在Data文件夹）
% 动力学问题：
%   载荷的时间函数可在SRC\Solvers\GeneralizedAlphaSolver.m中修改，默认是阶跃输入
%   也可改为正弦输入等

    % 设置路径 
    addpath(genpath('SRC')); 
    addpath(genpath('Data')); 
    
    % 输出文件，.out文件仅用于简单查看平动位移，.plt文件用于绘图
    [pathstr, name, ~] = fileparts(inputFileName);
    outputFileName = fullfile(pathstr, [name, '.out']);
    % pltFileName = fullfile(pathstr, [name, '.plt']);
    
    fidOut = fopen(outputFileName, 'w');
    if fidOut == -1, error('无法创建输出文件: %s', outputFileName); end
    finishup = onCleanup(@() fclose(fidOut));

    % 读取数据 
    femDomain = Domain.Instance();
    if ~femDomain.ReadData(inputFileName)
        error('读取输入文件失败。');
    end

    fprintf(fidOut, ' C O N T R O L   I N F O R M A T I O N\n\n');
    fprintf(fidOut, '      NUMBER OF NODAL POINTS . . . . . . . (NUMNP)  = %10d\n', femDomain.NUMNP);
    fprintf(fidOut, '      NUMBER OF ELEMENT GROUPS . . . . . . (NUMEG)  = %10d\n', femDomain.NUMEG);
    fprintf(fidOut, '      NUMBER OF LOAD CASES . . . . . . . . (NLCASE) = %10d\n', femDomain.NLCASE);
    fprintf(fidOut, '      SOLUTION MODE  . . . . . . . . . . . (MODEX)  = %10d\n\n', femDomain.MODEX);

    % 求解-
    tic;
    
    % 组装刚度矩阵
    femDomain.AssembleStiffnessMatrix();
    
    % 定义静力学和动力学结果的文件名，用于在Tecplot中画图
    fileStatic  = fullfile('Data', [name, '_Static.plt']);
    fileDynamic = fullfile('Data', [name, '_Dynamic.plt']);

    % 根据输入文件选择求解器
    if femDomain.MODEX == 1
    % 无论是否选择动力学，都先算一遍静力学作为“基准参考” 
        fprintf('\n [Step 1] Running Static Analysis for Reference...\n');
        staticSolver = StaticSolver();
        staticSolver.Solve(femDomain);
        fprintf('   Writing Static Reference to: %s\n', fileStatic);
        WriteTecplotLine(femDomain, fileStatic, 0.0, true);

        % 动力学求解
        if femDomain.AnalysisType == 1
            fprintf('\n [Step 2] Dynamic Analysis Selected. Starting...\n');
            p = femDomain.DynParams;
            solver = GeneralizedAlphaSolver(p.dt, p.nSteps, p.rho_inf, p.alpha, p.beta);
            solver.OutputFileName = fileDynamic;
            solver.Solve(femDomain);
            
        else
            fprintf('\n [Step 2] No Dynamic Analysis requested. Finished.\n');
        end
    end
    
        % solver.Solve(femDomain); 
    totalTime = toc;

    % 对于动力学，输出的最后一帧的位移；对于静力学，是平衡位置
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