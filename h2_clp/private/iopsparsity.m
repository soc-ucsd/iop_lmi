function const = iopsparsity(CYv,CUv,N,m,p,T_struct,R_struct)

%encode sparsity constraints

    const = [];
    z     = sym('z');


    %SPARSITY OF U
    Us = sym('U',[m p*(N+1)]);

    U_tf = zeros(m,p);
    for t = 1:N+1
        U_tf = U_tf + Us(:,[(t-1)*m+1:t*p])/z^(t-1);
    end

    for i=1:m      %ach1
        for j=1:p
            if(T_struct(i,j)==0)
                fprintf(' sparsity of U:  Percentage %6.4f \n', 100*(p*(i-1)+j)/m/p );
                [num,~] = numden(U_tf(i,j));
                cc      = coeffs(num,z);
                [A_eq,b_eq]    = equationsToMatrix(cc,vec(Us));
                A_eqs          = double(A_eq);
                b_eqs          = double(b_eq);
                const    = [const, A_eqs*vec(CUv)== b_eqs];
            end
        end
    end



    %SPARSITY OF Y
    Ys = sym('Y',[p p*(N+1)]);

    Y_tf = zeros(p,p);
    for t = 1:N+1
        Y_tf = Y_tf + Ys(:,[(t-1)*p+1:t*p])/z^(t-1);
    end



    for i=1:p      
        for j=1:p
            if(R_struct(i,j)==0)
                fprintf(' sparsity of Y:  Percentage %6.4f \n', 100*(p*(i-1)+j)/m/p );
                [num,~] = numden(Y_tf(i,j));
                cc      = coeffs(num,z);
                [A_eq,b_eq]    = equationsToMatrix(cc,vec(Ys));
                A_eqs          = double(A_eq);
                b_eqs          = double(b_eq);
                const    = [const, A_eqs*vec(CYv)== b_eqs];
            end
        end
    end

end