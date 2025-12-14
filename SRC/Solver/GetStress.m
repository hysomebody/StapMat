%* *****************************************************************
%* - Function of STAPMAT in Solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     Calculation of strain and stress                            *
%*                                                                 *
%* - Call procedures:                                              *
%*     SRC/Mechanics/Truss/TrussStress.m   - TrussStress()         *       
%*                                                                 *
%* - Called by :                                                   *
%*     ./Solver.m                                                  *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.22                *
%*                                                                 *
%* *****************************************************************

function GetStress(NUM)

% Different type of element
global cdata;
NUMEG = cdata.NUMEG;
IOUT = cdata.IOUT;

for N = 1:NUMEG
    NPAR1 = cdata.NPAR(1);
    if (NPAR1 == 1) TrussStress(NUM, N)
    else error(' *** ERROR *** No Such Element'); end 
end



end