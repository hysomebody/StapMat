function WriteTecplotVolume(domain, filename, currentTime, isNewFile)
% WRITETECPLOTVOLUME 输出四面体单元结果到 Tecplot (含总位移，修复行过长问题)
% 适配 GeneralizedAlphaSolver 的时间步输出接口

    NUMNP = domain.NUMNP;
    
    NUME = 0;
    for g = 1:domain.NUMEG
        NUME = NUME + length(domain.ElemGroups{g});
    end
    
    NodalData = zeros(8, NUMNP); 
    
    U_temp = zeros(domain.NEQ, 1); 
    
    for i = 1:NUMNP
        node = domain.NodeList(i);
        d_vec = node.Displacement(1:3); % 仅取平动 U, V, W
        
        NodalData(1:3, i) = node.XYZ;          
        NodalData(4:6, i) = d_vec;
        NodalData(7, i)   = norm(d_vec); % 总位移
        
        if isprop(node, 'Temperature')
            NodalData(8, i) = node.Temperature;
        else
            NodalData(8, i) = 0;
        end
        
        % 填充全局 U_temp
        for dof=1:6
            eq = node.BCode(dof);
            if eq > 0, U_temp(eq) = node.Displacement(dof); end
        end
    end
    
    % 收集单元应力
    ElemStress  = zeros(1, NUME);
    eCount = 0;
    for g = 1:domain.NUMEG
        group = domain.ElemGroups{g};
        for e = 1:length(group)
            elem = group(e);
            eCount = eCount + 1;
            
            % 计算应力
            res = elem.CalcStress(U_temp);
            ElemStress(eCount) = res.Stress; 
        end
    end
    
    if isNewFile
        fid = fopen(filename, 'w');
        fprintf(fid, 'TITLE = "FEM Volume Output"\n');
        fprintf(fid, 'VARIABLES = "X", "Y", "Z", "U", "V", "W", "TotalDisp", "Temp", "VonMises"\n');
    else
        fid = fopen(filename, 'a');
    end
    
    if fid == -1, error('无法打开文件: %s', filename); end
    
    fprintf(fid, 'ZONE T="Time=%.4f", STRANDID=1, SOLUTIONTIME=%.6e, ', currentTime, currentTime);
    fprintf(fid, 'NODES=%d, ELEMENTS=%d, ZONETYPE=FETETRAHEDRON, DATAPACKING=BLOCK, ', NUMNP, NUME);
    fprintf(fid, 'VARLOCATION=([9]=CELLCENTERED)\n');
    
    for varIdx = 1:8
        dataVec = NodalData(varIdx, :);
        WriteDataWithLineBreaks(fid, dataVec, 10);
    end
    
    WriteDataWithLineBreaks(fid, ElemStress, 10);
    
    for g = 1:domain.NUMEG
        group = domain.ElemGroups{g};
        for e = 1:length(group)
            elem = group(e);
            fprintf(fid, '%d %d %d %d\n', ...
                elem.Nodes(1).ID, elem.Nodes(2).ID, elem.Nodes(3).ID, elem.Nodes(4).ID);
        end
    end
    
    fclose(fid);
end

% -------------------------------------------------------------------------
% 辅助函数：带换行的写入
% -------------------------------------------------------------------------
function WriteDataWithLineBreaks(fid, dataVec, numPerLine)
    n = length(dataVec);
    fmtChunk = repmat('%15.6e ', 1, numPerLine);
    fmtChunk = [fmtChunk, '\n']; 
    
    nChunks = floor(n / numPerLine);
    remainder = mod(n, numPerLine);
    
    if nChunks > 0
        blockData = reshape(dataVec(1:nChunks*numPerLine), numPerLine, nChunks);
        fprintf(fid, fmtChunk, blockData); 
    end
    
    if remainder > 0
        idxStart = nChunks * numPerLine + 1;
        for i = idxStart:n
            fprintf(fid, '%15.6e ', dataVec(i));
        end
        fprintf(fid, '\n'); 
    end
end