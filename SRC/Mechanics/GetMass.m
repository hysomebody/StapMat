% GETMASS 组装结构整体一致质量矩阵 (Consistent Mass Matrix)
% 用法:
%   直接调用 GetMass()。
%   成功执行后，结果存储在全局变量 sdata.M_Global 中 (稀疏矩阵格式)。
%
% 注意:
%   必须确保单元数据已读取 (sdata.MATP, sdata.E, sdata.RHO 等不为空)。
%   建议在 ReadFile 阶段完成数据读取。
%
% 适用场景:
%   1. 动力学分析 (SolveDynamic)
%   2. 特征值分析 (计算自振频率)
%   3. 静力学自重计算 (Gravity Load)
%
% -Programmed by:
%   宋博
% -----------------------------------------------------------
function GetMass()
    global cdata sdata;
    % --- 防御性检查 ---
    if ~isfield(sdata, 'RHO') && cdata.NPAR(1) == 2 % 如果是梁单元
        error('错误: 未找到密度数据 sdata.RHO。请检查 ReadBeam 是否已在 ReadFile 中运行。');
    end
    fprintf('Forming Mass Matrix...\n');

    NEQ = sdata.NEQ;
    NUMEG = cdata.NUMEG;

    % 为了高效组装稀疏矩阵，先收集三元组 (i, j, value)
    % 预分配内存，这里预估每行大约有20个非零元
    estimated_nz = NEQ * 20; 
    I_idx = zeros(estimated_nz, 1);
    J_idx = zeros(estimated_nz, 1);
    V_val = zeros(estimated_nz, 1);
    count = 0;

    % 循环遍历所有单元组
    for Group = 1:NUMEG
        % 获取该组单元的类型
        Type = cdata.NPAR(1); % 这里假设只有一个组或简单的NPAR结构，STAPMAT原版NPAR可能在循环里读
        % 注意：STAPMAT原版处理多组单元的逻辑在 GetStiff 里是读文件循环。
        % 为了简化，我们假设 sdata 已经存储了所有单元信息。
        % 如果你的 ReadBeam 把所有单元都存在 sdata.MATP 里了，我们直接遍历单元即可。
        
        % 获取该组的单元总数 (这里简化逻辑，假设 sdata 中存储了所有单元列表)
        % 实际上 STAPMAT 的存储结构比较散，我们需要根据 ReadBeam 的存储方式来。
        % 基于之前的 ReadBeam，sdata.MATP 存储了所有单元的材料ID，sdata.LM 存储了自由度。
        
        NUME = size(sdata.MATP, 1); % 单元总数
        
        for N = 1:NUME
            MTYPE = sdata.MATP(N);
            % 判断是否为梁单元 (假设梁单元的特征是 NDOF=6 或根据 Type 判断)
            % 这里简单处理：只要是 ReadBeam 读进去的，且自由度为12，就按梁处理
            
            % 获取材料属性
            E = sdata.E(MTYPE);
            Nu = sdata.Nu(MTYPE);
            A = sdata.AREA(MTYPE);
            Iy = sdata.Iy(MTYPE);
            Iz = sdata.Iz(MTYPE);
            J_val = sdata.J_T(MTYPE);
            RHO = sdata.RHO(MTYPE); % 密度
            
            % 获取几何信息 (同 BeamStiff)
            XI = sdata.XYZ(1:3, N);
            XJ = sdata.XYZ(4:6, N);
            K_ID = sdata.K_Node(N);
            XK = [sdata.X(K_ID); sdata.Y(K_ID); sdata.Z(K_ID)];
            
            % 计算坐标变换
            v1 = XJ - XI; 
            L = norm(v1);
            if L < 1e-10, continue; end
            v1 = v1 / L;
            v_ref = XK - XI;
            v3 = cross(v1, v_ref); v3 = v3/norm(v3);
            v2 = cross(v3, v1); v2 = v2/norm(v2);
            Lam = [v1'; v2'; v3'];
            T = blkdiag(Lam, Lam, Lam, Lam);
            
            % --- 核心：调用单元质量矩阵函数 ---
            % 注意：必须确保 BeamMass.m 在路径中
            m_elem = BeamMass(RHO, A, Iy, Iz, J_val, L, T);
            
            % --- 收集组装索引 ---
            lm_vec = sdata.LM(:, N);
            
            for i = 1:12
                row = lm_vec(i);
                if row > 0
                    for j = 1:12
                        col = lm_vec(j);
                        if col > 0
                            count = count + 1;
                            % 动态扩容 (如果预估小了)
                            if count > length(I_idx)
                                I_idx(end+1 : end+NEQ*10) = 0;
                                J_idx(end+1 : end+NEQ*10) = 0;
                                V_val(end+1 : end+NEQ*10) = 0;
                            end
                            
                            I_idx(count) = row;
                            J_idx(count) = col;
                            V_val(count) = m_elem(i, j);
                        end
                    end
                end
            end
        end
    end

    % 截断多余的零
    I_idx = I_idx(1:count);
    J_idx = J_idx(1:count);
    V_val = V_val(1:count);

    % 生成全局稀疏质量矩阵
    sdata.M_Global = sparse(I_idx, J_idx, V_val, NEQ, NEQ);
    
    fprintf('Mass Matrix Assembled. Size: %d x %d\n', NEQ, NEQ);
end