% FILE PATH: SRC/Solvers/GeneralizedAlphaSolver.m
classdef GeneralizedAlphaSolver < Solver
    % GENERALIZEDALPHASOLVER 3rd-Order Generalized-Alpha Method
    % Reference: Dynamics.txt (genalpha3)
    
    properties
        dt          % Time step
        NSteps      % Number of time steps
        rho_inf     % Spectral radius
        RayleighA   % Rayleigh Damping alpha
        RayleighB   % Rayleigh Damping beta
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
            
            % 3. Algorithm Parameters [cite: 3-5]
            rho = obj.rho_inf;
            alpha_m = (13 + 20*rho - 5*rho^2) / (12 * (rho + 1)^2);
            alpha_f = (1 + 3*rho) / (2 * (rho + 1)^2);
            gamma1  = 5/12 + alpha_m - alpha_f;
            gamma2  = gamma1; 
            
            % 4. System Matrices (E, F) [cite: 5-6]
            I = speye(nd);
            Z = sparse(nd, nd);
            E_mat = [I, Z; Z, M];
            F_mat = [Z, -I; K, C];
            
            % 5. Effective Stiffness [cite: 15]
            dt_val = obj.dt;
            Keff = (dt_val * alpha_m) * E_mat + (dt_val^2 * alpha_f * gamma1) * F_mat;
            
            fprintf('   Factorizing Effective Stiffness Matrix...\n');
            try
                KeffDec = decomposition(Keff, 'lu');
            catch
                [L_fac, U_fac, P_fac, Q_fac] = lu(Keff);
            end
            
            % =============================================================
            % [修改重点] 定义载荷函数句柄 fhandle，整合协作者的逻辑
            % =============================================================
            
            % A. 获取空间分布 (从输入文件读取 Load Case 1)
            F_static = domainObj.AssembleForce(1);
            
            % B. 定义时间函数 (平滑加载 Ramp, 1.0秒达到峰值)
            %    注意：为了防止 t=0 时导数爆炸，必须使用 Ramp 而非阶跃
            time_func = @(t) min(t / 1.0, 1.0) .* (t >= 0); 
            % 如果你想用正弦波验证 (算例2)，可以取消注释下面这行：
            % time_func = @(t) sin(1.0 * t); 
            
            % C. 组合成 fhandle (这是协作者代码的核心接口)
            fhandle = @(t) F_static * time_func(t);
            
            % D. 定义 g_of_time 
            %    完全复刻 Dynamics.txt 的形式
            function gt = g_of_time(tt)
                ft = fhandle(tt);
                gt = [zeros(nd, 1); ft];
            end
            
            % =============================================================
            
            % 6. Initial Conditions
            u0 = zeros(nd, 1);
            v0 = zeros(nd, 1);
            
            % Initial Acceleration a0 [cite: 9]
            f0 = fhandle(0); % 调用 fhandle
            a0 = M \ (f0 - C*v0 - K*u0);
            
            x     = [u0; v0];
            x_dot = [v0; a0];
            
            % Initialize x_ddot (Using Finite Difference) [cite: 10-12]
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
                
                % Calculate RHS using g_of_time [cite: 17-18]
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
                
                % Update State [cite: 20-21]
                V1 = V_vec + dt_val*(1 - gamma1)*A_vec + dt_val*gamma1*A1;
                U1 = U_vec + dt_val*V_vec + (dt_val^2/2)*(1 - gamma2)*A_vec + (dt_val^2/2)*gamma2*A1;
                
                x      = U1;
                x_dot  = V1;
                x_ddot = A1;
                
                if mod(n, 100) == 0
                    fprintf('Step %d/%d done.\n', n, obj.NSteps);
                end
            end
            
            % 8. Post-Processing
            u_final = x(1:nd); 
            domainObj.UpdateNodalDisplacements(u_final);
            
            fprintf(' [GeneralizedAlphaSolver] Analysis Complete.\n');
        end
    end
end