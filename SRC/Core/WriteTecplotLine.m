function WriteTecplotLine(domain, filename, currentTime, isNewFile)
% WRITETECPLOTLINE 输出梁/杆单元结果到 Tecplot 
%

    NUMNP = domain.NUMNP;
    
    % 统计总单元数
    NUME = 0;
    for g = 1:domain.NUMEG
        NUME = NUME + length(domain.ElemGroups{g});
    end
    
    % 所有节点的坐标和位移 [X, Y, Z, U, V, W, Rx, Ry, Rz]
    NodalData = zeros(9, NUMNP); 

    NUME = 0;
    for g = 1:domain.NUMEG
        NUME = NUME + length(domain.ElemGroups{g});
    end

    % 临时的全局位移向量，计算应力用
    U_temp = zeros(domain.NEQ, 1); 
    
    for i = 1:NUMNP
        node = domain.NodeList(i);
        NodalData(1:3, i) = node.XYZ;          
        NodalData(4:9, i) = node.Displacement; 
        
        % 填充 U_temp
        for dof=1:6
            eq = node.BCode(dof);
            if eq > 0, U_temp(eq) = node.Displacement(dof); end
        end
    end
    
    % 单元数据：Stress, Moment, Torsion
    ElemStress  = zeros(1, NUME);
    ElemMoment  = zeros(1, NUME); 
    ElemTorsion = zeros(1, NUME); 
    eCount = 0;
    for g = 1:domain.NUMEG
        group = domain.ElemGroups{g};
        for e = 1:length(group)
            elem = group(e);
            eCount = eCount + 1;

            res = elem.CalcStress(U_temp);
            ElemStress(eCount) = res.Stress; 
            ElemMoment(eCount)  = res.Moment; 
            ElemTorsion(eCount) = res.Torsion; 
        end
    end
    
    % 写入文件
    if isNewFile
        fid = fopen(filename, 'w');
        % 文件头
        fprintf(fid, 'TITLE = "FEM Line Output"\n');
        fprintf(fid, 'VARIABLES = "X", "Y", "Z", "U", "V", "W", "RX", "RY", "RZ", "Stress", "Moment", "Torsion"\n');
    else
        fid = fopen(filename, 'a');
    end
    
    if fid == -1, error('无法打开文件: %s', filename); end
    
    % Zone 头
    fprintf(fid, 'ZONE T="Time=%.4f", STRANDID=1, SOLUTIONTIME=%.6e, ', currentTime, currentTime);
    fprintf(fid, 'NODES=%d, ELEMENTS=%d, ZONETYPE=FELINESEG, DATAPACKING=BLOCK, ', NUMNP, NUME);
    fprintf(fid, 'VARLOCATION=([10,11,12]=CELLCENTERED)\n');
    
    % 数据块（1-9）
    for varIdx = 1:9
        fprintf(fid, '%15.6e ', NodalData(varIdx, :));
        fprintf(fid, '\n');
    end
    
    % 单元变量 (10-12)
    fprintf(fid, '%15.6e ', ElemStress);  fprintf(fid, '\n'); % Stress
    fprintf(fid, '%15.6e ', ElemMoment);  fprintf(fid, '\n'); % Moment
    fprintf(fid, '%15.6e ', ElemTorsion); fprintf(fid, '\n'); % Torsion
    
    % 连接关系
    for g = 1:domain.NUMEG
        group = domain.ElemGroups{g};
        for e = 1:length(group)
            elem = group(e);
            n1 = elem.Nodes(1).ID;
            n2 = elem.Nodes(2).ID;
            fprintf(fid, '%d %d\n', n1, n2);
        end
    end
    
    fclose(fid);
end