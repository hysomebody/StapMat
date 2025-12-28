function ReadTe4()

    global sdata;
    sdata.NNODE = 4;
    sdata.NDOF = 3;

    % Read Material information
    ReadMaterial()

    % Read Element information
    ReadElements()

    % the second time stamp
    global cdata;
    cdata.CTIM(2,:) = clock;

end

% ----------------------- Functions -----------------------------------
% Read Material information
function ReadMaterial()

    global cdata;
    global sdata;
    % Get file pointers
    IIN = cdata.IIN;
    IOUT = cdata.IOUT;

    if (sdata.NPAR(3) == 0) sdata.NPAR(3) = 1; end
    fprintf(IOUT, '\n M A T E R I A L   D E F I N I T I O N\n');
    fprintf(IOUT, '\n NUMBER OF DIFFERENT SETS OF MATERIAL\n');
    fprintf(IOUT, ' AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . = %10d\n', ...
        sdata.NPAR(3));
    fprintf(IOUT, '  SET       YOUNG''S        Possion\n');
    fprintf(IOUT, ' NUMBER     MODULUS          ratio\n');
    fprintf(IOUT, '               E              NU\n');


    % Read material datas
    sdata.NUMMAT = sdata.NPAR(3);
    NUMMAT = sdata.NPAR(3);
    sdata.RHO = zeros(NUMMAT, 1);
    sdata.E = zeros(NUMMAT, 1);
    sdata.NU = zeros(NUMMAT, 1);
    sdata.D = zeros(NUMMAT, 6 , 6);
    sdata.AREA = ones(NUMMAT, 1);
    sdata.ALPHA = ones(NUMMAT, 1);
    for I = 1:sdata.NPAR(3)
        tmp = str2num(fgetl(IIN));
        N = round(tmp(1));
        sdata.RHO(N) = tmp(2);
        sdata.E(N) = tmp(3);
        sdata.NU(N) = tmp(4);
        sdata.ALPHA = tmp(5);
        A1 = tmp(4)/(1-tmp(4));
        A2 = (1-2*tmp(4))/2/(1-tmp(4));
        A3 = tmp(3)*(1-tmp(4))/(1+tmp(4))/(1-2*tmp(4));
        D = [ 1  A1 A1 0  0  0
              A1  1 A1 0  0  0
              A1 A1  1 0  0  0
               0  0  0 A2 0  0
               0  0  0 0  A2 0
               0  0  0 0  0 A2];
        D = D*A3;
        sdata.D(N,:,:) = D;
        fprintf(IOUT, '%5d    %12.5e  %14.6e\n', N, tmp(3), tmp(4));
    end


end

% Read elements information
function ReadElements()

    global cdata;
    global sdata;

    % Get file pointer
    IIN = cdata.IIN;
    IOUT = cdata.IOUT;

    fprintf(IOUT, '\n\n E L E M E N T   I N F O R M A T I O N\n');
    fprintf(IOUT, '\n      ELEMENT          NODE          NODE          NODE          NODE       MATERIAL\n');
    fprintf(IOUT, '      NUMBER-N           I             J             K             L         SET NUMBER\n');

    % Get Position data
    NUME = sdata.NPAR(2);NLCASE = sdata.NLCASE;
    XYZ = zeros(sdata.NNODE*sdata.NDOF, NUME);
    MATP = zeros(1, NUME);                                     % the type of material
    LM = zeros(sdata.NNODE*sdata.NDOF, NUME);                  % connectivity matrix
    X = sdata.X; Y = sdata.Y; Z = sdata.Z; ID = sdata.ID;TE =sdata.TE;TN=sdata.TN;
    XYZ = sdata.XYZ; 
    for N = 1:NUME
        tmp = str2num(fgetl(IIN));
        I = round(tmp(2));
        J = round(tmp(3));
        K = round(tmp(4));
        L = round(tmp(5));
        MTYPE = round(tmp(6));

    %   Save element information
        XYZ(1, N) = X(I);
        XYZ(2, N) = Y(I);
        XYZ(3, N) = Z(I);
        XYZ(4, N) = X(J);
        XYZ(5, N) = Y(J);
        XYZ(6, N) = Z(J);
        XYZ(7, N) = X(K);
        XYZ(8, N) = Y(K);
        XYZ(9, N) = Z(K);
        XYZ(10, N) = X(L);
        XYZ(11, N) = Y(L);
        XYZ(12, N) = Z(L);
        TE(1,N)=TN(I);
        TE(2,N)=TN(J);
        TE(3,N)=TN(K);
        TE(4,N)=TN(L);
        MATP(N) = MTYPE;

        fprintf(IOUT, '%10d      %10d    %10d    %10d    %10d       %5d\n', N, I, J, K, L, MTYPE);

    %   Compute connectivity matrix
        LM(1, N) = ID(1, I);
        LM(4, N) = ID(1, J);
        LM(2, N) = ID(2, I);
        LM(5, N) = ID(2, J);
        LM(3, N) = ID(3, I);
        LM(6, N) = ID(3, J);
        LM(7, N) = ID(1, K);
        LM(10, N) = ID(1, L);
        LM(8, N) = ID(2, K);
        LM(11, N) = ID(2, L);
        LM(9, N) = ID(3, K);
        LM(12, N) = ID(3, L);

    end
    sdata.XYZ = [sdata.XYZ,XYZ]; 
    sdata.LM = [sdata.LM,LM];
    sdata.MATP = [sdata.MATP,MATP];
    sdata.TE = [sdata.TE,TE];
    
    sdata.STRAIN = [sdata.STRAIN;zeros(NUME, NLCASE, 'double')];
    sdata.STRESS = [sdata.STRESS;zeros(NUME, NLCASE, 'double')];

    % Clear the memory of X, Y, Z
    sdata.X = double(0);
    sdata.Y = double(0);
    sdata.Z = double(0);


end