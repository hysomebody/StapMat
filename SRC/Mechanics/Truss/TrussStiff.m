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
sdata.NDOF = 3;

end

% Assemble structure stiffness matrix
function Assemble()
global sdata;
global cdata;
S = zeros(6, 6, 'double');
ST = zeros(6, 1, 'double');
sdata.STIFF = zeros(sdata.NWK, 1, 'double');

NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
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
    ST(4) = -ST(1); ST(5) = -ST(2); ST(6) = -ST(3);
    
    for J = 1:6
        YY = ST(J) * XX;
        for I = 1:J 
            S(I, J) = ST(I)*YY; 
        end
    end
    
%   SRC/Mechanics/ADDBAN.m
    ADDBAN(S, LM(:, N));
    
end

% The third time stamp
cdata.TIM(3, :) = clock;

end
