function BeamStress(NUM_LoadCase)
    global cdata sdata;
    IOUT = cdata.IOUT;
    NUME = cdata.NPAR(2);
    
    fprintf(IOUT, '\n\n BEAM FORCES AND MOMENTS (Load Case %d)\n', NUM_LoadCase);
    fprintf(IOUT, ' ELEM      Axial        Shear-Y        Shear-Z       Torque       Moment-Y       Moment-Z\n');
    
    for N = 1:NUME
        MTYPE = sdata.MATP(N);
        E = sdata.E(MTYPE);
        Nu = sdata.Nu(MTYPE);
        A = sdata.AREA(MTYPE);
        Iy = sdata.Iy(MTYPE);
        Iz = sdata.Iz(MTYPE);
        J_val = sdata.J_T(MTYPE);
        G = E / (2 * (1 + Nu));
        
        % 获取几何
        XI = sdata.XYZ(1:3, N);
        XJ = sdata.XYZ(4:6, N);
        K_ID = sdata.K_Node(N);
        XK = [sdata.X(K_ID); sdata.Y(K_ID); sdata.Z(K_ID)];
        
        v1 = XJ - XI; 
        L = norm(v1);
        v1 = v1 / L;
        v_temp = XK - XI;
        v3 = cross(v1, v_temp);
        v3 = v3 / norm(v3);
        v2 = cross(v3, v1);
        
        Lam = [v1'; v2'; v3']; 
        T = blkdiag(Lam, Lam, Lam, Lam);
        
        % 获取单元节点位移
        U_elem_global = zeros(12, 1);
        lm_vec = sdata.LM(:, N);
        
        % 从全局位移向量 DIS 中提取
        GlobalDis = sdata.DIS(:, NUM_LoadCase);
        for i = 1:12
            if lm_vec(i) > 0
                U_elem_global(i) = GlobalDis(lm_vec(i));
            end
        end
        
        % 转换到局部坐标系
        u_local = T * U_elem_global;
        
        % 计算内力 F = k_local * u_local
        k_local = BeamLocalStiff(E, G, A, Iy, Iz, J_val, L);
        F_local = k_local * u_local;
        
        % 提取节点 J 的内力作为单元内力输出 (或者输出两端)
        % F_local 顺序: [Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2]
        % 这里输出节点 J 端 (7-12) 的内力，符号通常取反
        Fx = F_local(7);
        Fy = F_local(8);
        Fz = F_local(9);
        Mx = F_local(10);
        My = F_local(11);
        Mz = F_local(12);
        
        fprintf(IOUT, '%5d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n', ...
            N, Fx, Fy, Fz, Mx, My, Mz);
    end
end