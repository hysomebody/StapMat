%* *****************************************************************
%* - Function of STAPMAT in initialization phase                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     Read input file of STAPMAT                                  *
%*                                                                 *
%* - Call procedures:                                              *
%*     SRC/Initiation/ReadFile.m - InitBasicData()                 *
%*                                                                 *
%* - Called by :                                                   *
%*     stapmat.m                                                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.21                *
%*                                                                 *
%* *****************************************************************
function ReadFile2(fname)
    fname = strcat('.\Data\', fname);  % Deal the filename
    
    % Get global class
    global cdata;
    global sdata;
    
    % Open files
    cdata.IIN = fopen(fname, 'r');
    
    % Begin Read input file
    fprintf('Input phase ...\n\n');
    
    % the first time stamp
    cdata.TIM = zeros(5, 6, 'double');
    cdata.TIM(1,:) = clock;
    
    IIN = cdata.IIN;
 
    %% Read Control data
    cdata.HED = fgetl(IIN);
    
    tmp = str2num(fgetl(IIN));
    cdata.NUMNP = int64(tmp(1));
    cdata.NUMEG = int64(tmp(2));
    cdata.NLCASE = int64(tmp(3));
    cdata.MODEX = int64(tmp(4));
    
    if (cdata.NUMNP == 0)
        return;
    end
    
    %% Read nodal point data
    % 检查 NDOF 是否已初始化 (防止忘记调用 InitBeam)
    if isempty(sdata) || ~isfield(sdata, 'NDOF') || sdata.NDOF == 0
         error('Error: sdata.NDOF is not set. Please call InitBeam() (or InitTruss) before ReadFile().');
    end

    InitBasicData();
    
    % Define local variables to speed up
    ID = sdata.ID; 
    X = sdata.X; 
    Y = sdata.Y; 
    Z = sdata.Z;
    
    for i = 1:cdata.NUMNP
        tmp = str2num(fgetl(IIN));
        
        % 动态读取 NDOF 个边界条件
        for j = 1:sdata.NDOF
            % tmp(1)是节点号，所以约束码从 tmp(2) 开始
            ID(j, i) = int64(tmp(1 + j));
        end
    
        % 坐标数据紧跟在边界条件之后
        % 索引位置 = 1(节点号) + NDOF(约束码数量) + 1(X坐标) 
        idx_coord = 1 + sdata.NDOF; 
        X(i) = double(tmp(idx_coord + 1));
        Y(i) = double(tmp(idx_coord + 2));
        Z(i) = double(tmp(idx_coord + 3));
    end
    
    sdata.ID = ID; 
    sdata.X = X; 
    sdata.Y = Y; 
    sdata.Z = Z;
    
    %% Compute the number of equations
    sdata.IDOrigin = ID;
    NEQ = 0;
    
    for N = 1:cdata.NUMNP
        for I = 1:sdata.NDOF  % 修改前写死了I = 1:3
            if (ID(I,N) == 0)
                NEQ = NEQ + 1;
                ID(I,N) = NEQ;
            else
                ID(I,N) = 0;
            end
        end
    end
    
    sdata.ID = ID;
    sdata.NEQ = NEQ;
    
    %% Read load data
    % Init control data
    NLCASE = cdata.NLCASE;
    sdata.R = zeros(NEQ, NLCASE, 'double');
    R = sdata.R;
    
    % Read data Loop
    for N = 1:cdata.NLCASE
        tmp = str2num(fgetl(IIN));
        cdata.LL = int64(tmp(1)); 
        cdata.NLOAD = int64(tmp(2));
        NLOAD = cdata.NLOAD;
        
        % Init load data
        sdata.NOD = zeros(NLOAD, 1, 'int64');
        sdata.IDIRN = zeros(NLOAD, 1, 'int64');
        sdata.FLOAD = zeros(NLOAD, 1, 'double');
        
        NOD = sdata.NOD; 
        IDIRN = sdata.IDIRN; 
        FLOAD = sdata.FLOAD;
        
        % Read load data lines
        for I = 1:NLOAD
            tmp = str2num(fgetl(IIN));
            NOD(I) = int64(tmp(1));
            IDIRN(I) = int64(tmp(2));
            FLOAD(I) = double(tmp(3));
        end
        
        if (cdata.MODEX == 0)
            return;
        end
        
        % Compute load vector
        for L = 1:NLOAD
            % 安全检查：防止 IDIRN 超过 NDOF
            if IDIRN(L) > sdata.NDOF
                error('Load direction (%d) exceeds NDOF (%d).', IDIRN(L), sdata.NDOF);
            end
            
            II = ID(IDIRN(L), NOD(L));
            if (II > 0)
                R(II, N) = R(II, N) + FLOAD(L);
            end
        end
        
        sdata.NOD = NOD; 
        sdata.IDIRN = IDIRN; 
        sdata.FLOAD = FLOAD; 
        sdata.R = R;
    end 
    % --- 荷载读取循环结束 ---
    
    %% Read Element Data 
    
    fprintf('Reading Element Data...\n');
    
    % 记录当前文件位置，为了处理 Truss 的兼容性问题
    current_fpos = ftell(IIN);
    
    for N = 1:cdata.NUMEG
        % 读取这一组单元的控制行 (例如: 2 2 1 -> 类型2, 2个单元, 1种材料)
        line_str = fgetl(IIN);
        tmp = str2num(line_str);
        
        % 存入 cdata.NPAR 供后续使用
        cdata.NPAR = tmp; 
        ElementType = cdata.NPAR(1);
        
        if ElementType == 1
            % --- TRUSS 处理逻辑 ---
            fprintf('  Found Truss Group (Type 1). Skipping read in ReadFile to preserve compatibility.\n');
            
            % 关键操作：因为我们不想改 TrussStiff，所以这里不能把数据读走。
            % 我们把文件指针退回到读这一行之前的位置！
            % 这样等会儿调用 GetStiff -> TrussStiff 时，它能读到这一行。
            fseek(IIN, current_fpos, 'bof'); 
            
            % 注意：如果这里回退了，for 循环里的下一次 fgetl 还会读到这一行，导致死循环。
            % 所以如果决定不读 Truss，这里必须 break 或者 return，或者你要非常清楚后续逻辑。
            % 鉴于 STAPMAT 结构，Element 数据通常在文件最后。
            % 如果这里不读，就让 ReadFile 结束即可。
            return; 
            
        elseif ElementType == 2
            % --- BEAM 处理逻辑 ---
            fprintf('  Reading Beam Group (Type 2)...\n');
            ReadBeam(); % 调用 ReadBeam.m，它会把该组所有数据读完
            
            % ReadBeam 读完后，更新文件位置标记，准备读下一组（如果有）
            current_fpos = ftell(IIN);
        end
    end

end

%% Functions

% InitBasicData
function InitBasicData()
    global cdata;
    global sdata;
    
    cdata.NPAR = zeros(10, 1, 'int64');
    
    % 使用 sdata.NDOF 动态分配内存
    sdata.ID = zeros(sdata.NDOF, cdata.NUMNP, 'int64'); 
    sdata.X = zeros(cdata.NUMNP, 1, 'double');
    sdata.Y = zeros(cdata.NUMNP, 1, 'double');
    sdata.Z = zeros(cdata.NUMNP, 1, 'double');
end