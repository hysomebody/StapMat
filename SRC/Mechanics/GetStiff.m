%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Forming the stiffness matrix                                *
%*                                                                 *
%* - Call procedures:                                              *
%*     SRC/Mechanics/Truss/TrussStiff.m - TrussStiff()             *
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

function GetStiff()
% Get global variables
global cdata;

% Read the type of element
IIN = cdata.IIN;
IOUT = cdata.IOUT;
fprintf(IOUT, '\n\n E L E M E N T   G R O U P   D A T A\n');

for N = 1:cdata.NUMEG
    if (N ~= 1) fprintf(IOUT,'\n'); end
    tmp = str2num(fgetl(IIN));
    for I = 1:length(tmp) cdata.NPAR(I) = tmp(I); end

    fprintf(IOUT, '\n\n E L E M E N T   D E F I N I T I O N\n');
    fprintf(IOUT, ['\n ELEMENT TYPE  . . . . . . . . . . . . .( NPAR(1) ) . . = %10d\n' ...
        '     EQ.1, TRUSS ELEMENTS\n' ...
        '     EQ.2, ELEMENTS CURRENTLY\n' ...
        '     EQ.3, NOT AVAILABLE\n' ...
        ' NUMBER OF ELEMENTS. . . . . . . . . . .( NPAR(2) ) . . = %10d\n'], ...
        cdata.NPAR(1), cdata.NPAR(2));

%   Different kinds of element
    NPAR1 = cdata.NPAR(1);
    if (NPAR1 == 1) TrussStiff()
    else error(' *** ERROR *** No Such Element'); end
    
end

end
