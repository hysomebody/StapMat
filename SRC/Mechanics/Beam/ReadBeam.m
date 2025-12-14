function ReadBeam()
    % 读取材料和单元信息
    ReadBeamMaterial();
    ReadBeamElements();
    
    global cdata;
    cdata.TIM(2,:) = clock; % 记录时间
end

function ReadBeamMaterial()
    global cdata sdata;
    IIN = cdata.IIN;
    IOUT = cdata.IOUT;
    
    fprintf(IOUT, '\n BEAM MATERIAL DEFINITION \n');
    fprintf(IOUT, ' MAT_ID       E           NU          AREA          Iy            Iz            J\n         Density');
    
    NUMMAT = cdata.NPAR(3);
    % 初始化数组
    sdata.E = zeros(NUMMAT, 1);
    sdata.Nu = zeros(NUMMAT, 1);
    sdata.AREA = zeros(NUMMAT, 1);
    sdata.Iy = zeros(NUMMAT, 1);
    sdata.Iz = zeros(NUMMAT, 1);
    sdata.J_T = zeros(NUMMAT, 1);
    sdata.RHO = zeros(NUMMAT, 1);
    
    for i = 1:NUMMAT
        line = fgetl(IIN);
        tmp = str2num(line);
        if isempty(tmp), continue; end
        
        N = round(tmp(1));
        sdata.E(N) = tmp(2);
        sdata.Nu(N) = tmp(3);
        sdata.AREA(N) = tmp(4);
        sdata.Iy(N) = tmp(5);
        sdata.Iz(N) = tmp(6);
        sdata.J_T(N) = tmp(7);
        sdata.RHO(N) = tmp(8);

        fprintf(IOUT, '%5d %12.5e %10.3f %12.5e %12.5e %12.5e %12.5e %12.5e\n', ...
            N, sdata.E(N), sdata.Nu(N), sdata.AREA(N), sdata.Iy(N), sdata.Iz(N), sdata.J_T(N), sdata.RHO(N));
    end
end

function ReadBeamElements()
    global cdata sdata;
    IIN = cdata.IIN;
    IOUT = cdata.IOUT;
    
    NUME = cdata.NPAR(2);
    % 初始化数据
    sdata.XYZ = zeros(3*2, NUME); % 2个节点的坐标
    sdata.MATP = zeros(NUME, 1);
    sdata.K_Node = zeros(NUME, 1); % 存储方向节点K
    sdata.LM = zeros(12, NUME);    % 12个自由度 (2节点 * 6 DOF)
    
    X = sdata.X; Y = sdata.Y; Z = sdata.Z; ID = sdata.ID;
    
    fprintf(IOUT, '\n BEAM ELEMENT DEFINITION \n');
    fprintf(IOUT, ' ELEM     NODE_I    NODE_J    NODE_K    MAT_ID\n');
    
    for N = 1:NUME
        line = fgetl(IIN);
        tmp = str2num(line);
        if isempty(tmp), continue; end
        
        Num = round(tmp(1));
        I = round(tmp(2));
        J = round(tmp(3));
        K = round(tmp(4)); % 方向节点
        MTYPE = round(tmp(5));
        
        % 存储几何信息
        sdata.XYZ(1,N)=X(I); sdata.XYZ(2,N)=Y(I); sdata.XYZ(3,N)=Z(I);
        sdata.XYZ(4,N)=X(J); sdata.XYZ(5,N)=Y(J); sdata.XYZ(6,N)=Z(J);
        
        sdata.MATP(N) = MTYPE;
        sdata.K_Node(N) = K;
        
        % 生成LM数组 (Connectivity Matrix)
        % 节点I的 6 个DOF
        for k = 1:6
            sdata.LM(k, N) = ID(k, I);
        end
        % 节点J的 6 个DOF
        for k = 1:6
            sdata.LM(k+6, N) = ID(k, J);
        end
        
        fprintf(IOUT, '%5d %8d %8d %8d %8d\n', Num, I, J, K, MTYPE);
    end
    
    % 更新带宽等辅助数据
    ColHt(sdata.LM);  % 需要确保 ColHt 支持 12 行的 LM
end