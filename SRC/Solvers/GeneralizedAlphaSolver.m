classdef GeneralizedAlphaSolver < Solver
    % 广义-alpha 法动力学求解器 (3阶精度)
    % by 聂学峰
    
    properties
        dt          % 时间步长
        NSteps      % 总时间步数
        rho_inf     % 谱半径
        RayleighA   % 阻尼系数 α
        RayleighB   % 瑞利阻尼系数 β

        % 结果存储
        % TimeHistory % 存储时间点
        % DispHistory % 存储特定节点的位移历史 (用于绘图对比)
        % StressHistory % 存储特定单元的应力历史
        OutputFileName % 存储 Tecplot 输出文件名
    end
    
    methods
        function obj = GeneralizedAlphaSolver(dt, N, rho, a, b)
            obj.dt = dt;
            obj.NSteps = N;
            obj.rho_inf = rho;
            if nargin < 4, a = 0.0; end
            if nargin < 5, b = 0.001; end
            obj.RayleighA = a;
            obj.RayleighB = b;
        end
        
        function Solve(obj, domainObj)
            fprintf(' [GeneralizedAlphaSolver] Start Dynamic Analysis...\n');
            
            if isempty(domainObj.GlobalM), domainObj.AssembleMassMatrix(); end
            M = domainObj.GlobalM;
            K = domainObj.GlobalK;
            nd = domainObj.NEQ;

            % 阻尼矩阵 C = alpha*M + beta*K
            C = obj.RayleighA * M + obj.RayleighB * K;
            fprintf('   Damping Matrix Constructed (Rayleigh: a=%.4f, b=%.4f)\n', obj.RayleighA, obj.RayleighB);
            
            rho = obj.rho_inf;
            alpha_m = (13 + 20*rho - 5*rho^2) / (12 * (rho + 1)^2);
            alpha_f = (1 + 3*rho) / (2 * (rho + 1)^2);
            gamma1  = 5/12 + alpha_m - alpha_f;
            gamma2  = gamma1; 
            
            I = speye(nd);
            Z = sparse(nd, nd);
            E_mat = [I, Z; Z, M];% 根据文献构造
            F_mat = [Z, -I; K, C]; % 根据文献构造
             
            dt_val = obj.dt;
            Keff = (dt_val * alpha_m) * E_mat + (dt_val^2 * alpha_f * gamma1) * F_mat;
            
            fprintf('   Factorizing Effective Stiffness Matrix...\n');
            try
                KeffDec = decomposition(Keff, 'lu');
            catch
                [L_fac, U_fac, P_fac, Q_fac] = lu(Keff);
            end
            
            % 
            F_static = domainObj.AssembleForce(1);

            %--------------------------------------------------------
            % 定义时间函数 （输入的载荷在这里修改）
            % % 正弦载荷
            % omega = 60; 
            % time_func = @(t) sin(omega * t);

            % 阶跃载荷
            time_func = @(t) 1.0;
            
            % 组合成总载荷句柄 f(t)
            fhandle = @(t) F_static * time_func(t);

           
          % % A. 尝试监控第一个受载节点
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
          %  % C. 监控单元
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

            % 扩充后的载荷向量函数 g(t) 
            function gt = g_of_time(tt)
                ft = fhandle(tt);
                gt = [zeros(nd, 1); ft];
            end
            
            u0 = zeros(nd, 1);
            v0 = zeros(nd, 1);
            
            f0 = fhandle(0); 
            a0 = M \ (f0 - C*v0 - K*u0);
            
            x     = [u0; v0];
            x_dot = [v0; a0];

            epsFD = 1e-6;
            g0 = g_of_time(0);
            g_eps = g_of_time(epsFD);
            g_dot0 = (g_eps - g0) / epsFD;
            
            rhs0 = g_dot0 - F_mat * x_dot;
            x_ddot = E_mat \ rhs0;
            
            % 时间积分循环
            fprintf('   Time Integration: %d steps, dt=%.2e\n', obj.NSteps, dt_val);
            % 在命令行窗口显示进度
            
            for n = 1:obj.NSteps
                t_n = (n-1)*dt_val;
                t_naf = t_n + alpha_f * dt_val; 
                
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
                if ~isempty(obj.OutputFileName) && (mod(n, 10) == 0) % 每10步输出一帧，避免文件过大
                     % 必须先更新 Domain 中的位移，WriteTecplotLine 才能读到正确数据
                     domainObj.UpdateNodalDisplacements(x(1:nd));
                     
                     % 如果是第10步(第一次输出)，设为新文件(true)，否则为追加(false)
                     isFirstWrite = (n == 10); 
                     WriteTecplotLine(domainObj, obj.OutputFileName, t_naf, isFirstWrite);
                end

                % %  输出某个节点的位移（或其他量），对比用
                % csvFileName = 'History_Node82.csv';
                % 
                % if mod(n, 10) == 0
                % 
                %     targetNode = 8;
                %     dofPerNode = 6; 
                %     targetDir  = 3; 
                %     globalDofIndex = (targetNode - 1) * dofPerNode + targetDir;
                %     uv_val = x(globalDofIndex);
                % 
                %     if n == 10 
                %         % 第一次输出
                %         fid = fopen(csvFileName, 'w');
                %         if fid == -1, error('无法打开 CSV 文件写入'); end
                %         fprintf(fid, 'Time, Disp_Y_Node82\n');
                %     else
                %         % 后续输出
                %         fid = fopen(csvFileName, 'a');
                %         if fid == -1, error('无法打开 CSV 文件追加'); end
                %     end
                % 
                %     % 写入数据
                %     fprintf(fid, '%.6f, %.6e\n', t_naf, uv_val);
                %     fclose(fid);
                % end
                % 
                % if mod(n, 10) == 0
                %     fprintf('Step %d/%d done. Max Disp=%.2e\n', n, obj.NSteps, max(abs(x(1:nd))));
                % end
            end
            
            % % 简单的绘图验证 
            % figure(1);
            % subplot(2,1,1); plot(obj.TimeHistory, obj.DispHistory); 
            % title('Node 51 Displacement (STAPMAT)'); grid on;
            % subplot(2,1,2); plot(obj.TimeHistory, obj.StressHistory);
            % title('Element 50 Stress (STAPMAT)'); grid on;
            
            u_final = x(1:nd); 
            domainObj.UpdateNodalDisplacements(u_final);
            
            fprintf(' [GeneralizedAlphaSolver] Analysis Complete.\n');
        end
    end
end