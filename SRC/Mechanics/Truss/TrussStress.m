%* *****************************************************************
%* - Function of STAPMAT in solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     To calculate stresses                                       *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Solver/GetStress.m                                      *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.22                *
%*                                                                 *
%* *****************************************************************

function TrussStress(NUM, NG)

% Get global data
global cdata;
global sdata;

IOUT = cdata.IOUT;
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ;
E = sdata.E; AREA = sdata.AREA; LM = sdata.LM;
U = sdata.DIS(:, NUM);

fprintf(IOUT, ['\n\n  S T R E S S  C A L C U L A T I O N S  F O R  ' ...
    'E L E M E N T  G R O U P %4d\n\n' ... 
    '       ELEMENT             FORCE            STRESS\n' ...
    '       NUMBER\n'], NG);

for N = 1:NUME
    MTYPE = MATP(N);
   
%   compute the length of truss element
    DX = XYZ(1, N) - XYZ(4, N);
    DY = XYZ(2, N) - XYZ(5, N);
    DZ = XYZ(3, N) - XYZ(6, N);
    XL2 = DX*DX + DY*DY + DZ*DZ;
    
    ST(1) = DX / XL2 * E(MTYPE);
    ST(2) = DY / XL2 * E(MTYPE);
    ST(3) = DZ / XL2 * E(MTYPE);
    ST(4) = -ST(1); ST(5) = -ST(2); ST(6) = -ST(3);
    
    STR = 0.0;
    
    if (LM(1, N) > 0) STR = STR + ST(1)*U(LM(1, N)); end
    if (LM(2, N) > 0) STR = STR + ST(2)*U(LM(2, N)); end
    if (LM(3, N) > 0) STR = STR + ST(3)*U(LM(3, N)); end
    if (LM(4, N) > 0) STR = STR + ST(4)*U(LM(4, N)); end
    if (LM(5, N) > 0) STR = STR + ST(5)*U(LM(5, N)); end
    if (LM(6, N) > 0) STR = STR + ST(6)*U(LM(6, N)); end
    
    P = STR*AREA(MTYPE);
    
    fprintf(IOUT, ' %10d           %13.6e    %13.6e\n', N, P, STR);
end

end