% FILE PATH: SRC/Solvers/GeneralizedAlphaSolver.m
classdef GeneralizedAlphaSolver < Solver
    % GENERALIZEDALPHASOLVER 3rd-Order Generalized-Alpha Method
    % Reference: Dynamics.txt (genalpha3)
    
    properties
        dt          % Time step(时间步长)
        NSteps      % Number of time steps(步数)
        rho_inf     % Spectral radius 
        RayleighA   % Rayleigh Damping alpha (瑞利商 α)
        RayleighB   % Rayleigh Damping beta (瑞利商 β)

        % 结果存储
        % TimeHistory % 存储时间点
        % DispHistory % 存储特定节点的位移历史 (用于绘图对比)
        % StressHistory % 存储特定单元的应力历史
        OutputFileName % 存储 Tecplot 输出文件名
    end
    
    methods
        % Constructor
        function obj = GeneralizedAlphaSolver(dt, N, rho, a, b)
            obj.dt = dt;
            obj.NSteps = N;
            obj.rho_inf = rho;
            if nargin < 4, a = 0.0; end
            if nargin < 5, b = 0.001; end
            obj.RayleighA = a;
            obj.RayleighB = b;
        end
        
        % Execute Solution
        function Solve(obj, domainObj)
            fprintf(' [GeneralizedAlphaSolver] Start Dynamic Analysis...\n');
            
            % 1. Prepare Matrices
            if isempty(domainObj.GlobalM), domainObj.AssembleMassMatrix(); end
            M = domainObj.GlobalM;
            K = domainObj.GlobalK;
            nd = domainObj.NEQ;
            
            % 2. Damping Matrix
            C = obj.RayleighA * M + obj.RayleighB * K;
            fprintf('   Damping Matrix Constructed (Rayleigh: a=%.4f, b=%.4f)\n', obj.RayleighA, obj.RayleighB);
            
            % 3. Algorithm Parameters 
            rho = obj.rho_inf;
            alpha_m = (13 + 20*rho - 5*rho^2) / (12 * (rho + 1)^2);
            alpha_f = (1 + 3*rho) / (2 * (rho + 1)^2);
            gamma1  = 5/12 + alpha_m - alpha_f;
            gamma2  = gamma1; 
            
            % 4. System Matrices (E, F) 
            I = speye(nd);
            Z = sparse(nd, nd);
            E_mat = [I, Z; Z, M];
            F_mat = [Z, -I; K, C];
            
            % 5. Effective Stiffness 
            dt_val = obj.dt;
            Keff = (dt_val * alpha_m) * E_mat + (dt_val^2 * alpha_f * gamma1) * F_mat;
            
            fprintf('   Factorizing Effective Stiffness Matrix...\n');
            try
                KeffDec = decomposition(Keff, 'lu');
            catch
                [L_fac, U_fac, P_fac, Q_fac] = lu(Keff);
            end
            
            % 1. 获取空间分布 (从输入文件读取 Load Case 1)
            F_static = domainObj.AssembleForce(1);

            % 2. 定义时间函数 (针对算例: sin(omega * t))
            omega = 60; 
            time_func = @(t) sin(omega * t);
            
            % 3. 组合成 fhandle 
            fhandle = @(t) F_static * time_func(t);
            
          %  % 初始化历史记录
          %  % A. 尝试监控第一个受载节点
          %  lc = domainObj.LoadCases(1);
          %  if ~isempty(lc.Nodes)
          %      monitorNodeID = lc.Nodes(1);
          %  else
          %      % 无载荷时监控中间节点
          %      monitorNodeID = round(domainObj.NUMNP / 2);
          %  end
          %  
          %  % B. 边界检查
          %  if monitorNodeID > domainObj.NUMNP
          %       monitorNodeID = 1; 
          %  end
          %  % C. 监控单元 (取中间号)
          %  totalElements = 0;
          %  for g = 1:domainObj.NUMEG
          %      totalElements = totalElements + length(domainObj.ElemGroups{g});
          %  end
          %  monitorElemID = round(totalElements / 2);
          %  if monitorElemID > totalElements, monitorElemID = 1; end
          %  
          %  fprintf('   [Monitor] Tracking Node %d and Element %d\n', monitorNodeID, monitorElemID);
          %  
          %  obj.TimeHistory = zeros(obj.NSteps, 1);
          %  obj.DispHistory = zeros(obj.NSteps, 1); 
          %  obj.StressHistory = zeros(obj.NSteps, 1);

            % 4. 定义 g_of_time 
            function gt = g_of_time(tt)
                ft = fhandle(tt);
                gt = [zeros(nd, 1); ft];
            end
            
            % =============================================================
            
            % 5.初始化位移
            u0 = zeros(nd, 1);
            v0 = zeros(nd, 1);
            
            % Initial Acceleration a0 
            f0 = fhandle(0); % 调用 fhandle
            a0 = M \ (f0 - C*v0 - K*u0);
            
            x     = [u0; v0];
            x_dot = [v0; a0];
            
            % Initialize x_ddot (Using Finite Difference) 
            epsFD = 1e-6;
            g0 = g_of_time(0);
            g_eps = g_of_time(epsFD);
            g_dot0 = (g_eps - g0) / epsFD;
            
            rhs0 = g_dot0 - F_mat * x_dot;
            x_ddot = E_mat \ rhs0;
            
            % 7. Time Integration Loop
            fprintf('   Time Integration: %d steps, dt=%.2e\n', obj.NSteps, dt_val);
            
            for n = 1:obj.NSteps
                t_n = (n-1)*dt_val;
                t_naf = t_n + alpha_f * dt_val; % [cite: 16]
                
                U_vec = x;
                V_vec = x_dot;
                A_vec = x_ddot;
                
                % Calculate RHS using g_of_time 
                g_naf = g_of_time(t_naf);
                
                term1 = E_mat * V_vec;
                term2 = (dt_val * (1 - alpha_m)) * (E_mat * A_vec);
                term3 = F_mat * U_vec;
                term4 = dt_val * (F_mat * V_vec);
                term5 = (dt_val^2 * alpha_f * (1 - gamma1)) * (F_mat * A_vec);
                
                ERHS = g_naf - term1 - term2 - term3 - term4 - term5;
                
                if exist('KeffDec', 'var')
                    A1 = KeffDec \ ERHS;
                else
                    A1 = Keff \ ERHS;
                end
                
                % Update State 
                V1 = V_vec + dt_val*(1 - gamma1)*A_vec + dt_val*gamma1*A1;
                U1 = U_vec + dt_val*V_vec + (dt_val^2/2)*(1 - gamma2)*A_vec + (dt_val^2/2)*gamma2*A1;
                
                x      = U1;
                x_dot  = V1;
                x_ddot = A1;
                
                % 输出到 Tecplot
                % 只有当定义了文件名，且满足输出间隔时才写入
                if ~isempty(obj.OutputFileName) && (mod(n, 100) == 0) % 每100步输出一帧，避免文件过大
                     % 必须先更新 Domain 中的位移，WriteTecplotLine 才能读到正确数据
                     domainObj.UpdateNodalDisplacements(x(1:nd));
                     
                     % 如果是第100步(第一次输出)，设为新文件(true)，否则为追加(false)
                     isFirstWrite = (n == 100); 
                     WriteTecplotLine(domainObj, obj.OutputFileName, t_naf, isFirstWrite);
                end

               % % --- 记录数据 ---
               % if n <= obj.NSteps
               %     obj.TimeHistory(n) = t_naf; % 或者 t_n，取决于想记录哪个时刻
               %     
               %     % 记录监控节点的位移 (假设 X 方向 DOF=1)
               %     % 获取 Node 51 的全局方程号
               %     mNode = domainObj.NodeList(monitorNodeID);
               %     eqNum = mNode.BCode(1);
               %     if eqNum > 0
               %         obj.DispHistory(n) = x(eqNum);
               %     end
               %     
               %     % --- 计算并记录应力 (Task 3) ---
               %     % 为了性能，不计算所有单元，只计算监控单元
               %     % 如果需要全场应力，建议只在特定步数计算
               %     elem = domainObj.ElemGroups{1}(monitorElemID); % 假设 Group 1 是 Truss
               %     
               %     % 提取当前时刻的位移向量 U_vec (即 x(1:nd)) 传给单元
               %     % TrussElement.CalcStress 需要全局 U
               %     stressRes = elem.CalcStress(x(1:nd));
               %     obj.StressHistory(n) = stressRes.Stress;
               % end

                if mod(n, 100) == 0
                    fprintf('Step %d/%d done. Max Disp=%.2e\n', n, obj.NSteps, max(abs(x(1:nd))));
                end
            end
            
            % % 简单的绘图验证 
            % figure(1);
            % subplot(2,1,1); plot(obj.TimeHistory, obj.DispHistory); 
            % title('Node 51 Displacement (STAPMAT)'); grid on;
            % subplot(2,1,2); plot(obj.TimeHistory, obj.StressHistory);
            % title('Element 50 Stress (STAPMAT)'); grid on;
            
            % 8. Post-Processing
            u_final = x(1:nd); 
            domainObj.UpdateNodalDisplacements(u_final);
            
            fprintf(' [GeneralizedAlphaSolver] Analysis Complete.\n');
        end
    end
end