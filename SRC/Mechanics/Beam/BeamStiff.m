function BeamStiff()
    % 1. 初始化与读取
    InitBeam(); 
    ReadBeam();
    
    fprintf('Forming Stiffness Matrix for Beams...\n');
    
    % 2. 计算带宽 (Address)
    Addres(); 
    
    global cdata sdata;
    if cdata.MODEX == 0   
        return;
    end % Data Check 模式直接返回
    
    % 3. 组装刚度矩阵
    AssembleBeam();
    
    cdata.TIM(3,:) = clock;
end

function AssembleBeam()
    global sdata;
    NUME = cdata.NPAR(2);
    
    % 循环所有单元
    for N = 1:NUME
        MTYPE = sdata.MATP(N);
        
        % 材料属性
        E = sdata.E(MTYPE);
        Nu = sdata.Nu(MTYPE);
        A = sdata.AREA(MTYPE);
        Iy = sdata.Iy(MTYPE);
        Iz = sdata.Iz(MTYPE);
        J_val = sdata.J_T(MTYPE);
        G = E / (2 * (1 + Nu));
        
        % 节点 I, J 坐标
        XI = sdata.XYZ(1:3, N);
        XJ = sdata.XYZ(4:6, N);
        
        % 方向节点 K 坐标
        K_ID = sdata.K_Node(N);
        % 这里需要从原始全局坐标取 K 点，或者在 ReadBeam 时存下来
        XK = [sdata.X(K_ID); sdata.Y(K_ID); sdata.Z(K_ID)];
        
        % 计算局部坐标轴向量
        v1 = XJ - XI; 
        L = norm(v1);
        if L < 1e-10, error('Element %d length is zero!', N); end
        v1 = v1 / L;
        
        v_temp = XK - XI;
        v3 = cross(v1, v_temp); % 局部 z 轴 (垂直于 I-J-K 平面)
        if norm(v3) < 1e-6, error('Node I, J, K are collinear!'); end
        v3 = v3 / norm(v3);
        
        v2 = cross(v3, v1); % 局部 y 轴 (在平面内)
        v2 = v2 / norm(v2);

        % 构造变换矩阵 T (3x3 -> 12x12)
        Lam = [v1'; v2'; v3']; % 3x3 方向余弦
        T = blkdiag(Lam, Lam, Lam, Lam); % 12x12
        
        % 计算局部刚度矩阵 k (12x12)
        k_local = BeamLocalStiff(E, G, A, Iy, Iz, J_val, L);
        
        % 转换到全局坐标系
        k_global = T' * k_local * T;
        
        % 组装到整体刚度矩阵
        ADDBAN(k_global, sdata.LM(:, N));
    end
end