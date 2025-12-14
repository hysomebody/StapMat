%* *****************************************************************
%* - Function of STAPMAT in Solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     Perform L*D*L(T) factorization of stiffness matrix          *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     ./Solve.m                                                   *
%*                                                                 *
%* - Programmed in Fortran 90 by Xiong Zhang                       *
%*                                                                 *
%* - Adapted to Matlab by:                                         *
%*     Yan Liu, Computational Dynamics Group, School of Aerospace  *
%*     Engineering, Tsinghua University, 2019.02.22                *
%*                                                                 *
%* *****************************************************************

function LDLTFactor()

% Get global data
global sdata;

A = sdata.STIFF; 
MAXA = sdata.MAXA; NEQ = sdata.NEQ; 

for N = 1:NEQ
    KN = MAXA(N);
    KL = KN + 1;
    KU = MAXA(N + 1) - 1;
    KH = KU - KL;
    
    if (KH > 0)
        K = N - KH;
        IC = 0;
        KLT = KU;
        for J = 1:KH
            IC = IC + 1;
            KLT = KLT - 1;
            KI = MAXA(K);
            ND = MAXA(K+1) - KI - 1;
            if (ND > 0)
                KK = min(IC, ND);
                C = 0.0;
                for L = 1:KK C = C+A(KI+L)*A(KLT+L); end
                A(KLT) = A(KLT) - C;
            end
            K = K + 1;
        end
    end
    
    if (KH >= 0)
        K = N;
        B = 0.0;
        for KK = KL:KU
            K = K - 1;
            KI = MAXA(K);
            C = A(KK) / A(KI);
            B = B + C*A(KK);
            A(KK) = C;
        end
        A(KN) = A(KN) - B;
    end
    
    if (A(KN) <= 0)
        error(['STOP - Stiffness matrix is not positive definite\n' ...
            'Nonpositive number for equation %8d is %20.12e\n'], N, A(KN));
    end
end

sdata.STIFF = A;

end