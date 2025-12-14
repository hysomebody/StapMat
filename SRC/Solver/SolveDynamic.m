function SolveDynamic()
% SOLVEDYNAMIC 使用 Newmark-beta 方法求解动力学方程
% M * u_dd + C * u_d + K * u = F(t)

    global cdata sdata;
    
    % === 1. 动力学参数设置 (由于输入文件没给，这里手动设置) ===
    dt = 0.01;          % 时间步长 (秒)
    T_total = 1.0;      % 总时间 (秒)
    steps = round(T_total / dt);
    
    % 阻尼参数 (Rayleigh Damping: C = alpha*M + beta*K)
    alpha = 0.0;        % 质量比例阻尼
    beta  = 0.001;      % 刚度比例阻尼 (通常设一个小值)
    
    % Newmark 方法参数 (平均加速度法: 无条件稳定)
    gamma = 0.5;
    delta = 0.25;       % 有些书上写作 beta，为避免混淆这里叫 delta
    
    fprintf('\n>>> Starting Dynamic Analysis (Newmark Method)...\n');
    fprintf('    Time step: %.4f s, Total steps: %d\n', dt, steps);

    % === 2. 准备矩阵 ===
    NEQ = sdata.NEQ;
    
    % [A] 获取质量矩阵 M (稀疏矩阵)
    GetMass(); 
    M = sdata.M_Global;
    
    % [B] 获取刚度矩阵 K (从 Skyline 转换为 稀疏矩阵)
    % STAPMAT 通常有 Stiff2Sparse 函数，如果没有，你需要确认 SPSTIFF 是否已生成
    % 如果 MODEX=2 (MODEX在输入文件定义)，STAPMAT可能已经生成了 SPSTIFF
    % 这里我们强制转换一次以防万一
    if isfield(sdata, 'SPSTIFF') && ~isempty(sdata.SPSTIFF)
        K = sdata.SPSTIFF;
    else
        % 调用 STAPMAT 自带的转换工具 (假设叫 Stiff2Sparse)
        % 如果没有，你需要利用 IIndex, JIndex 和 STIFF 数组自己组装
        try
            K = Stiff2Sparse(); 
        catch
            error('Cannot convert Stiffness Matrix to Sparse. Check Stiff2Sparse function.');
        end
    end
    
    % [C] 构造阻尼矩阵 C
    C = alpha * M + beta * K;
    
    % === 3. 初始化状态向量 ===
    u = zeros(NEQ, 1);      % 位移
    v = zeros(NEQ, 1);      % 速度
    a = zeros(NEQ, 1);      % 加速度
    
    % 存储历史结果用于输出
    % DIS_History: 行=自由度, 列=时间步
    sdata.DIS_History = zeros(NEQ, steps);
    sdata.Time_History = zeros(1, steps);
    
    % === 4. 计算有效刚度矩阵 K_hat ===
    % Newmark 公式: K_hat = K + a0 * M + a1 * C
    a0 = 1 / (delta * dt^2);
    a1 = gamma / (delta * dt);
    a2 = 1 / (delta * dt);
    a3 = 1 / (2 * delta) - 1;
    a4 = gamma / delta - 1;
    a5 = dt / 2 * (gamma / delta - 2);
    a6 = dt * (1 - gamma);
    a7 = gamma * dt;
    
    K_hat = K + a0 * M + a1 * C;
    
    % 对 K_hat 进行分解 (预处理，提高后续求解速度)
    [L_fac, U_fac, P_fac] = lu(K_hat); % LU 分解
    
    % === 5. 初始加速度计算 ===
    % F0 = R (假设初始时刻就有力) - C*v0 - K*u0
    % 简单起见，假设初始静止 u=0, v=0
    % Load Vector R (静态荷载)
    R_static = sdata.R(:, 1); % 取第一个工况的荷载
    
    % 假设荷载是阶跃函数 (Step Load): t>0 时 F = R_static
    F_current = R_static; 
    
    % a0 = M \ (F0 - C*v0 - K*u0); 
    % 注意：这里用反斜杠求解
    a = M \ (F_current - C*v - K*u);
    
    % === 6. 时间步进循环 ===
    for i = 1:steps
        t = i * dt;
        
        % 1. 计算有效荷载 R_hat
        % 假设荷载保持恒定 (Step Load)。如果是地震波，这里要随时间变。
        F_next = R_static; 
        
        R_hat = F_next + M * (a0 * u + a2 * v + a3 * a) ...
                       + C * (a1 * u + a4 * v + a5 * a);
        
        % 2. 求解位移 u_next
        % K_hat * u_next = R_hat
        u_next = U_fac \ (L_fac \ (P_fac * R_hat));
        
        % 3. 更新加速度和速度
        a_next = a0 * (u_next - u) - a2 * v - a3 * a;
        v_next = v + a6 * a + a7 * a_next;
        
        % 4. 记录结果
        sdata.DIS_History(:, i) = u_next;
        sdata.Time_History(i) = t;
        
        % 5. 更新状态
        u = u_next;
        v = v_next;
        a = a_next;
        
        if mod(i, 10) == 0
            fprintf('    Step %d / %d done.\n', i, steps);
        end
    end
    
    fprintf('Dynamic Analysis Completed.\n');
    
    % 可以在这里调用 WriteDisDynamic 或者输出到 Tecplot
end