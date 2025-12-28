function Te4Stress(NUM, NG)

    % Get global data
    global cdata;
    global sdata;

    IOUT = cdata.IOUT;
    E1 = sdata.NUME(NG)+1;
    E2 = sdata.NUME(NG+1);
    MATP = sdata.MATP; 
    XYZ = sdata.XYZ;
    E = sdata.E; NU = sdata.NU; LM = sdata.LM;D=sdata.D;TE=sdata.TE;ALPHA=sdata.ALPHA;
    U = sdata.DIS(:, NUM);
    S = zeros(E2-E1+1,1);

    for N = E1:E2
        MTYPE = MATP(N);

        x=[XYZ(1, N),XYZ(4, N),XYZ(7, N),XYZ(10, N)];
        y=[XYZ(2, N),XYZ(5, N),XYZ(8, N),XYZ(11, N)];
        z=[XYZ(3, N),XYZ(6, N),XYZ(9, N),XYZ(12, N)];
        
        V6 = det([1 x(1) y(1) z(1) ;
                  1 x(2) y(2) z(2) ;
                  1 x(3) y(3) z(3) ;
                  1 x(4) y(4) z(4)]);
        
        B=zeros(6,12);

        for i=1:4
            j = mod(i,4)+1;
            m = mod(i+1,4)+1;
            p = mod(i+2,4)+1;
            bi = -det( [ 1 y(j) z(j) 
                         1 y(m) z(m) 
                         1 y(p) z(p)] );
            ci =  det( [ 1 x(j) z(j) 
                         1 x(m) z(m) 
                         1 x(p) z(p)] );
            di = -det( [ 1 x(j) y(j) 
                         1 x(m) y(m) 
                         1 x(p) y(p)] );
            B(:,(i-1)*3+1:(i-1)*3+3)=(-1)^(i+1)*[bi  0  0
                                                  0 ci  0
                                                  0  0 di
                                                  0 di ci
                                                 di  0 bi
                                                 ci bi  0];
        end
        B=B/V6;

        UN=zeros(12,1);
        for i=1:12
            if (LM(i, N) > 0) UN(i)=U(LM(i, N)); end
        end

        SN=squeeze(D(MTYPE,:,:))*(B*UN-ALPHA*[1 1 1 0 0 0]'*(sum(TE(:,N)))/4);
        a=SN(1)+SN(2)+SN(3);
        b=SN(1)*SN(2)+SN(2)*SN(3)+SN(3)*SN(1)-SN(4)^2-SN(5)^2-SN(6)^2;
        c=SN(1)*SN(2)*SN(3)+2*SN(4)*SN(5)*SN(6)-(SN(1)*SN(4)^2+SN(2)*SN(5)^2+SN(3)*SN(6)^2);
        s=roots([1 -a b -c]);
        STR=sqrt(s(1)^2+s(2)^2+s(3)^2-s(1)*s(2)-s(2)*s(3)-s(3)*s(1));

        S(N-E1+1) = STR;        
        fprintf(IOUT, ' %10d           %13.6e\n', N, STR);
    end

    sdata.STRESS(E1:E2, NUM) = S;
end