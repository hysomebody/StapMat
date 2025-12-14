function InitBeam()
    global sdata;
    % 设置梁单元节点数为2
    sdata.NNODE = 2;
    % 关键：将全局自由度提升为 6
    sdata.NDOF = 6; 
end