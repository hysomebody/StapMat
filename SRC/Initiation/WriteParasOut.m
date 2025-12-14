%* *****************************************************************
%* - Function of STAPMAT in initialization phase                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     Write parameters to output file of STAPMAT                  *
%*                                                                 *
%* - Call procedures: None                                         *
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

function WriteParasOut()
global cdata;
global sdata;
% Open file
cdata.IOUT = fopen('.\Data\STAPMAT.OUT', 'w');
IOUT = cdata.IOUT;

fprintf(IOUT, ['\n %s \n\n'...
    ' C O N T R O L   I N F O R M A T I O N\n\n'...
    '      NUMBER OF NODAL POINTS . . . . . . . . . . (NUMNP)  = %10d \n' ...
    '      NUMBER OF ELEMENT GROUPS . . . . . . . . . (NUMEG)  = %10d \n' ...
    '      NUMBER OF LOAD CASES . . . . . . . . . . . (NLCASE) = %10d \n' ...
    '      SOLUTION MODE  . . . . . . . . . . . . . . (MODEX)  = %10d \n' ...
    '         EQ.0, DATA CHECK \n' ...
    '         EQ.1, EXECUTION'], ...
    cdata.HED, cdata.NUMNP, cdata.NUMEG, cdata.NLCASE, cdata.MODEX);

% Write complete nodal data
ID = sdata.IDOrigin; X = sdata.X; Y = sdata.Y; Z = sdata.Z;
fprintf(IOUT, '\n\n\n N O D A L   P O I N T   D A T A \n\n');
fprintf(IOUT, '       NODE          BOUNDARY                         NODAL POINT\n');
fprintf(IOUT, '      NUMBER     CONDITION  CODES                     COORDINATES\n');
fprintf(IOUT, '                    X    Y    Z               X            Y            Z\n');
for i = 1:cdata.NUMNP
    fprintf(IOUT, '%10d      %5d%5d%5d      %13.3f%13.3f%13.3f\n', ...
        i, ID(1,i), ID(2,i), ID(3,i), X(i), Y(i), Z(i));
end
sdata.IDOrigin = 0;                              % Delete old ID array

% Write equation numbers
ID = sdata.ID;
fprintf(IOUT, '\n\n EQUATION NUMBERS\n');
fprintf(IOUT, '\n        NODE         DEGREES OF FREEDOM\n');
fprintf(IOUT, '       NUMBER\n');
fprintf(IOUT, '          N                  X         Y         Z\n');
for N=1:cdata.NUMNP
    fprintf(IOUT, ' %10d         %10d%10d%10d\n', ...
        N, ID(1,N), ID(2,N), ID(3,N));
end

% Write the load vector
% Only the first load vector
fprintf(IOUT, '\n\n L O A D   C A S E   D A T A\n');
LL = cdata.LL; NLOAD = cdata.NLOAD;
NOD = sdata.NOD; IDIRN = sdata.IDIRN; FLOAD = sdata.FLOAD;
for I = 1:1 %cdata.NLCASE
    fprintf(IOUT, '\n     LOAD CASE NUMBER . . . . . . . = %10d\n', LL);
    fprintf(IOUT, '     NUMBER OF CONCENTRATED LOADS . = %10d\n', NLOAD);
    
    if (LL ~= I)
        error(' *** ERROR *** LOAD CASES ARE NOT IN ORDER');
    end
    fprintf(IOUT, '\n\n        NODE       DIRECTION      LOAD\n');
    fprintf(IOUT, '       NUMBER                   MAGNITUDE\n');
    
    for N = 1:NLOAD(I)
        fprintf(IOUT,'%10d         %4d       %12.5e',...
            NOD(N), IDIRN(N), FLOAD(N));
    end
end
end