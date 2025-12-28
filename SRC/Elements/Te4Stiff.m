function Te4Stiff(E1,E2)

    global cdata;
    global sdata;
    S = zeros(12, 12, 'double');

    MATP = sdata.MATP; XYZ = sdata.XYZ; 
    E = sdata.E; NU = sdata.NU; LM = sdata.LM;RHO = sdata.RHO;D = sdata.D;TE=sdata.TE;ALPHA=sdata.ALPHA;
    R=sdata.R;
    
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
            bi = -det( [ 1 y(j) z(j) ;
                         1 y(m) z(m) ;
                         1 y(p) z(p)] );
            ci =  det( [ 1 x(j) z(j) ;
                         1 x(m) z(m) ;
                         1 x(p) z(p)] );
            di = -det( [ 1 x(j) y(j) ;
                         1 x(m) y(m) ;
                         1 x(p) y(p)] );
            B(:,(i-1)*3+1:(i-1)*3+3)=(-1)^(i+1)*[bi  0  0
                                                  0 ci  0
                                                  0  0 di
                                                  0 di ci
                                                 di  0 bi
                                                 ci bi  0];
        end
        B=B/V6;

        V=abs(V6/6);
        
        f=ALPHA*B'*squeeze(D(MTYPE,:,:))*[1 1 1 0 0 0]'*(sum(TE(:,N)))/4*V;
        fe=[f;f;f;f];
        for NC = 1:sdata.NLCASE
            for i=1:12
                if (LM(i, N) > 0) R(LM(i, N),NC)=R(LM(i, N),NC)+fe(i); end
            end
        end
        sdata.R=R;
        
        S=B'*squeeze(D(MTYPE,:,:))*B*V;



        AddBoundary(S, LM(:, N), 1);

    end
end

