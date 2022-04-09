function Constraints = iopcons(G,CYv,CUv,CWv,CZv,N,z,opts)
% encoding achivability constraint

    [p,m] = size(G);      % system dimension
    CYs = sym('CX',[p p*(N+1)]);
    CUs = sym('CY',[m p*(N+1)]);
    CWs = sym('CW',[p m*(N+1)]);
    CZs = sym('CZ',[m m*(N+1)]);
    
    %%Express X,Y,W,Z as FIR transfer matrices of order N
    Y = zeros(p,p);
    U = zeros(m,p);
    W = zeros(p,m);
    Z = zeros(m,m);
    for t = 1:N+1
        Y = Y + CYs(:,[(t-1)*p+1:t*p])/z^(t-1);
        U = U + CUs(:,[(t-1)*p+1:t*p])/z^(t-1);
        W = W + CWs(:,[(t-1)*m+1:t*m])/z^(t-1);
        Z = Z + CZs(:,[(t-1)*m+1:t*m])/z^(t-1);
    end


    %%%achievability constraints
    ach1 = Y - G*U - eye(p);

    
    ach2 = W - G*Z;

    
    ach3 = -Y*G + W;

    ach4 = -U*G + Z - eye(m);
    
    
    Constraints=[];

    for i=1:p      %ach1
        for j=1:p
            
            fprintf(' ach1:  Percentage %6.4f \n', 100*(p*(i-1)+j)/p/p );
            [num,~] = numden(ach1(i,j));
            cc      = coeffs(num,z);
            [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CYs);vec(CUs)]);
            A_eqs          = double(A_eq);
            b_eqs          = double(b_eq);
            Constraints    = [Constraints, A_eqs*[vec(CYv);vec(CUv)]== b_eqs];
        end
    end

    for i = 1:p       %ach2
        for j = 1:m
            fprintf(' ach2:  Percentage %6.4f \n', 100*(m*(i-1)+j)/m/p );
            [num,~] = numden(ach2(i,j));
            cc      = coeffs(num,z);
            [A_eq,b_eq] = equationsToMatrix(cc,[vec(CWs);vec(CZs)]);
            A_eqs       = double(A_eq);
            b_eqs       = double(b_eq);
            Constraints = [Constraints, A_eqs*[vec(CWv);vec(CZv)]== b_eqs];
        end
    end

    for(i=1:p)       %ach3
        for(j=1:m)
            fprintf(' ach3:  Percentage %6.4f \n', 100*(m*(i-1)+j)/m/p );
            [num,~]=numden(ach3(i,j));
            cc=coeffs(num,z);
            [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CYs);vec(CWs)]);
            A_eqs   = double(A_eq);
            b_eqs    = double(b_eq);
            Constraints = [Constraints, A_eqs*[vec(CYv);vec(CWv)]== b_eqs];
        end
    end

    for(i=1:m)       %ach4
        for(j=1:m)
            fprintf(' ach4:  Percentage %6.4f \n', 100*(m*(i-1)+j)/m/m );
            [num,~]=numden(ach4(i,j));
            cc=coeffs(num,z);
            [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CUs);vec(CZs)]);
            A_eqs   = double(A_eq);
            b_eqs    = double(b_eq);
            Constraints = [Constraints, A_eqs*[vec(CUv);vec(CZv)]== b_eqs];
        end
    end
end
