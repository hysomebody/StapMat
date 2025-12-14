%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Compute the stiffness matrix of truss                       *
%*                                                                 *
%* - Call procedures:                                              *
%*     TrussStiff.m - InitTruss()                                  *
%*     ./ReadTruss.m - ReadTruss()                                 *
%*     SRC/Mechanics/Addres.m - Addres()                           *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Mechanics/GetStiff.m                                    *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.21                *
%*                                                                 *
%* *****************************************************************

function TrussStiff()

% Init variables of the element
InitTruss();

% Read Material and Elements
ReadTruss();

fprintf('Solution phase ...\n\n');

% calculate addresses of diagonal elements
Addres();

% Data check Or Solve
global cdata;
if (cdata.MODEX == 0) 
    cdata.TIM(3,:) = clock;
    cdata.TIM(4,:) = clock;
    cdata.TIM(5,:) = clock;
    return; 
end

% Assemble structure stiffness matrix
Assemble();




end

% ----------------------- Functions -----------------------------------

% Init parameters of truss element
function InitTruss()
global sdata;
sdata.NNODE = 2;
% sdata.NDOF = 3; % 全局自由度不能修改

end

% Assemble structure stiffness matrix

function Assemble()
global sdata;
global cdata;
% S = zeros(6, 6, 'double');
S_Expanded = zeros(12, 12, 'double'); %为了连接梁单元，桁架单元需要修改为12×12
ST = zeros(6, 1, 'double');
sdata.STIFF = zeros(sdata.NWK, 1, 'double');

%NUME = sdata.NUME;
NUME = cdata.NPAR(2);
MATP = sdata.MATP; XYZ = sdata.XYZ; 
E = sdata.E; AREA = sdata.AREA; LM = sdata.LM;
for N = 1:NUME
    MTYPE = MATP(N);
    
%   compute the length of truss element
    DX = XYZ(1, N) - XYZ(4, N);
    DY = XYZ(2, N) - XYZ(5, N);
    DZ = XYZ(3, N) - XYZ(6, N);
    XL2 = DX*DX + DY*DY + DZ*DZ;
    XL = sqrt(XL2);
    
    XX = E(MTYPE) * AREA(MTYPE) * XL;
    
    ST(1) = DX / XL2;
    ST(2) = DY / XL2;
    ST(3) = DZ / XL2;
    ST(4) = -ST(1); 
    ST(5) = -ST(2); 
    ST(6) = -ST(3);
    S_small = zeros(6, 6);
        for J = 1:6
            YY = ST(J) * XX;
            for I = 1:J
                S_small(I,J) = ST(I) * YY;
                S_small(J,I) = S_small(I,J); % 补全对称，方便下面映射
            end
        end
    
%   for J = 1:6
%       YY = ST(J) * XX;
%       for I = 1:J 
%           S(I, J) = ST(I)*YY; 
%       end
%   end

% 清空大矩阵
    S_Expanded(:) = 0;
    
% 定义映射关系：
% Truss的 1,2,3 (节点1平动) -> Beam体系的 1,2,3
% Truss的 4,5,6 (节点2平动) -> Beam体系的 7,8,9

% 块1: Node 1 平动 (对应 S_Expanded 1:3, 1:3)
    S_Expanded(1:3, 1:3) = S_small(1:3, 1:3);
    
% 块2: Node 1-2 耦合 (对应 S_Expanded 1:3, 7:9)
    S_Expanded(1:3, 7:9) = S_small(1:3, 4:6);
    
% 块3: Node 2-1 耦合 (对应 S_Expanded 7:9, 1:3)
    S_Expanded(7:9, 1:3) = S_small(4:6, 1:3);
    
% 块4: Node 2 平动 (对应 S_Expanded 7:9, 7:9)
    S_Expanded(7:9, 7:9) = S_small(4:6, 4:6);
    
% --- 修改点 3: 调用 ADDBAN ---
% 此时 sdata.LM(:, N) 应该是 12 行的 (因为 ReadFile 里根据 NDOF=6 生成了)
% S_Expanded 也是 12x12 的，完美匹配
    ADDBAN(S_Expanded, sdata.LM(:, N));

%%   SRC/Mechanics/ADDBAN.m
%    ADDBAN(S, LM(:, N));
    
end

% The third time stamp
cdata.TIM(3, :) = clock;

end
