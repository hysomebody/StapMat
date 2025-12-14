%* *****************************************************************
%* - Function of STAPMAT in Solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     To solve finite element static equilibrium equations        *
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

function ColSol(NUM)

% Get global data
global sdata;
A = sdata.STIFF; MAXA = sdata.MAXA; R = sdata.R(:,NUM);
NEQ = sdata.NEQ; NWK = sdata.NWK; 
NNM = NEQ + 1;

% Reduce right-hand-side load vector
for N = 1:NEQ
    KL = MAXA(N) + 1;
    KU = MAXA(N+1) - 1;
    if (KU-KL >= 0)
        K = N;
        C = 0.0;
        for KK = KL:KU
            K = K - 1;
            C = C + A(KK) * R(K);
        end
        R(N) = R(N) - C;
    end
end

% Back-Substitute
for N = 1:NEQ
    K = MAXA(N);
    R(N) = R(N) / A(K);
end

if (NEQ == 1) return; end;

N = NEQ;
for L = 2:NEQ
    KL = MAXA(N) + 1;
    KU = MAXA(N+1) - 1;
    if (KU-KL >= 0)
        K = N;
        for KK = KL:KU
            K = K - 1;
            R(K) = R(K) - A(KK)*R(N);
        end
    end
    N = N - 1;
end

sdata.DIS(:, NUM) = R(:);

end