function const = slssparsity(C,CNv,CLv,N,m,n,p,T_struct,R_struct)
%encode sparsity constraints for SLS


    const = [];
    CNs   = sym('CN',[n p*(N+1)]);
    z     = sym('z');


    N_tf = zeros(n,p);
    for t = 1:N+1
        N_tf = N_tf + CNs(:,[(t-1)*n+1:t*p])/z^(t-1);
    end

    %sparsity of C2N+I
    eq = eye(p)+C*N_tf;


    for i=1:p      %ach1
        for j=1:p
            if(R_struct(i,j)==0)
                fprintf(' sparsity of I+C2N:  Percentage %6.4f \n', 100*(p*(i-1)+j)/p/p );
                [num,~] = numden(eq(i,j));
                cc      = coeffs(num,z);
                [A_eq,b_eq]    = equationsToMatrix(cc,vec(CNs));
                A_eqs          = double(A_eq);
                b_eqs          = double(b_eq);
                const    = [const, A_eqs*vec(CNv)== b_eqs];
            end
        end
    end



    Ls = sym('L',[m p*(N+1)]);

    L_tf = zeros(m,p);
    for t = 1:N+1
        L_tf = L_tf + Ls(:,[(t-1)*p+1:t*p])/z^(t-1);
    end



    for i=1:m      
        for j=1:p
            if(T_struct(i,j)==0)
                fprintf(' sparsity of L:  Percentage %6.4f \n', 100*(m*(i-1)+j)/m/p );
                [num,~] = numden(Ls(i,j));
                cc      = coeffs(num,z);
                [A_eq,b_eq]    = equationsToMatrix(cc,vec(Ls));
                A_eqs          = double(A_eq);
                b_eqs          = double(b_eq);
                const    = [const, A_eqs*vec(CLv)== b_eqs];
            end
        end
    end
end