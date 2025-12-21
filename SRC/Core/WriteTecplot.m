%* *****************************************************************
%* - Function of STAPMAT in initialization phase                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     输出到Tecplot，包括节点位置、位移和网格                      *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     stapmat.m                                                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     cyz                                                         *
%*                                                                 *
%* *****************************************************************

function WriteTecplot(femDomain)

NLCASE = femDomain.NLCASE;
NUMNP = femDomain.NUMNP;
NUME = 0;

for L = 1:NLCASE

    timestamp = string(datetime("now",'Format', 'yyyyMMdd-HHmm'));
    fname = sprintf('%s-Ex%s-%s.plt', timestamp, string(femDomain.Title), string(L));
    fullpath = fullfile('.', 'Data', fname);
    
    % 确保 Data 目录存在
    if ~exist(fullfile('.', 'Data'), 'dir')
       mkdir(fullfile('.', 'Data'));
    end

    PLTOUT = fopen(fullpath, 'w');

    fprintf(PLTOUT,'TITLE = "Example : %s"\n',femDomain.Title);
    
    % 1. 修改变量头，增加 "Temp"
    fprintf(PLTOUT,'VARIABLES= "X"   "Y"   "Z"   "disp" "dispX"   "dispY"   "dispZ" "Temp"\n');

    % 总单元数
    for grp = 1:femDomain.NUMEG
        NUME = NUME + length(femDomain.ElemGroups{grp});
    end

    fprintf(PLTOUT,'Zone , N=%d, E=%d, ZONETYPE=FETETRAHEDRON, DATAPACKING=POINT\n', NUMNP, NUME);

    for II = 1:NUMNP
        node = femDomain.NodeList(II);
        d = node.Displacement;
        
        % 2. 输出节点数据时增加 node.Temperature
        % 格式：X_deformed, Y_deformed, Z_deformed, DispMag, u, v, w, Temp
        fprintf(PLTOUT,'%18.6e %18.6e %18.6e %18.6e %18.6e %18.6e %18.6e %18.6e\n', ...
            node.XYZ(1)+d(1), node.XYZ(2)+d(2), node.XYZ(3)+d(3), ...
            norm(d), d(1), d(2), d(3), node.Temperature);
    end

    % 遍历输出每个单元的节点编号
    for grp = 1:femDomain.NUMEG
        elements = femDomain.ElemGroups{grp};

        if isempty(elements)
            continue;
        end
        
        % 获取该组单元的节点数 (Beam/Truss=2, Tetra=4)
        NUM_Nodes = length(elements(1).Nodes);
        
        for e = 1:length(elements)
            el = elements(e);
            
            % Tecplot 的 FEQUADRILATERAL 需要 4 个节点。
            % 如果是 2 节点单元 (梁/杆)，重复最后节点以填充。
            
            for n=1:4
                % 简单的防止索引越界逻辑 (如果单元只有2个节点，ceil(3/2)=2, ceil(4/2)=2)
                idx = ceil(n * NUM_Nodes / 4.0); 
                if idx > NUM_Nodes, idx = NUM_Nodes; end
                
                fprintf(PLTOUT,'%18d', el.Nodes(idx).ID);
            end
            fprintf(PLTOUT,'\n');
        end

    end
    
    fclose(PLTOUT);
end



