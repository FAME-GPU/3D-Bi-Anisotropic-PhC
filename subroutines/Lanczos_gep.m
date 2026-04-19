function [EV,EW,nit] = Lanczos_gep(fncA,fncB,fncinvB,n,eigen_wanted,option)
    m         = option.m;     % size of expanded subspace
    k         = eigen_wanted;  % size of reduced subspace
    nit       = 0;
    itmax     = option.itmax;
    tolerence = option.tol;

    start_vec = option.v0;
    start_vec = start_vec/(option.v0'*fncB(option.v0));
    
    Q(:,1)    = start_vec;
    alpha     = zeros(m,1);
    beta      = zeros(m,1);
    [Q, alpha, beta] = Lanczos_step(fncA, fncB, fncinvB, n, Q, alpha, beta, 0, m);

    while nit < itmax
        nit = nit + 1;
        % expend Krylov subspace by Lanczos
        [Q, alpha, beta] = Lanczos_step(fncA, fncB, fncinvB, n, Q, alpha, beta, k, m);
        % restart step
        [Q, alpha, beta] = Implicit_restart_step(Q, alpha, beta, k, m);
        % stoping criterion
        while abs(beta(k)) < tolerence*abs(alpha(k))
            T       = diag(alpha(1:k)) + diag(beta(1:k-1), -1) + diag(beta(1:k-1), 1);
            [EV,EW] = eigen_solve(Q, T);
            return
        end
    end
end

%% Lanczos Step
function [Q, alpha, beta] = Lanczos_step(fncA, fncB, fncinvB, n, Q, alpha, beta, k, m)
    for j = k+1 : m
        temp     = fncA(Q(:,j));
        alpha(j) = real( Q(:,j)'*temp );
        if j == 1
            r = fncinvB(temp) - alpha(j)*Q(:,j);
        else
            r = fncinvB(temp) - alpha(j)*Q(:,j) - beta(j-1)*Q(:,j-1);
        end
        [Q(:,j+1), ~, beta(j)] = gsreorthog(Q(:,1:j), fncB , r, n, j);
    end
end
%% Implicit Restart step
function [Q, alpha, beta] = Implicit_restart_step(Q, alpha, beta, k, m)
    % construst matrix T
    T = diag(alpha(1:m)) + diag(beta(1:m-1),-1) + diag(beta(1:m-1),1);

    % choose unwanted ritz values
    unwanted_ritz = choose_unwanted(T,k,m);

    e       = zeros(m,1); 
    e(m)    = 1;
    I       = eye(m);
    U_total = I;
    
    % single shift step
    for i = 1:m-k
        [U,R]   = qr(T - unwanted_ritz(i)*I);
        U_total = U_total*U;
        T       = R*U + unwanted_ritz(i)*I;
        for j = 1:length(T)-2
            T(j,j+2:length(T)) = 0;
        end
        T = (T + T')/2;
    end
    b = U_total'*e;
    Q(:,1:m) = Q(:,1:m)*U_total;
    
    % new Arnoldi decomposition
    P          = T(k+1,k)*Q(:,k+1) + beta(m)*b(k)*Q(:,m+1);
    rho        = norm(P);
    Q          = Q(:,1:k);
    Q(:,k+1)   = P/rho;
    T          = T(1:k,1:k);
    beta(1:k)  = [diag(T,-1);rho];
    alpha(1:k) = diag(T);
end

function unwanted_ritz = choose_unwanted(T,k,m)
    ritz_val      = eig(T);
    [~,idx]       = sort(abs(ritz_val),'ascend');
    unwanted_ritz = ritz_val(idx(1:m-k));
end

%% Eigen Solver
function [EV,EW] = eigen_solve(Q,T)
    [EV, EW] = eig(T);
    EV       = Q(:,1:length(T))*EV; 
end
%% Reorthogonalize 
function [p, r, tho] = gsreorthog(Q, fncB, x, m ,n)
% [p, r, tho] = gsreorthog(Q, x)
%      x = Q*r + tho*p 
% Orthogonalize x against Q with reorthogonalization
    r = zeros(n,1);
    tho = x'*fncB(x);% tho = norm(x);
    if tho == 0
        p = zeros(n,1);
        return
    end
    mu = tho;
    x_bar = x;

    % find the index of the column with minimal 1-norm
    n1 = Q(:,1)'*fncB(Q(:,1));%n1 = norm(Q(:,1));
    k = 1;
    for i=2:n
        n2 = Q(:,i)'*fncB(Q(:,i));% n2 = norm(Q(:,i));
        if(n1>n2)
            n1 = n2; 
            k = i;
        end
    end

%     ei = zeros(m,1);
    ei = zeros(m,1);
    ei(k) = 1;
    temp = Q';
    while(1) 
        s = temp*fncB(x_bar);
%         s = Q'*x_bar;
        r = r + s;

        % x_bar is the residual
        x_bar = x_bar - Q*s; % x_bar = x_bar - Q*s;
        tau = x_bar'*fncB(x_bar);%tau = norm(x_bar);
        if(tau/tho >=1/2)
            break;
        end

        if(tau > 0.1*mu*eps)
            tho = tau;
        else
            tho = 0.1*tho*eps;
            mu = tho;
            x_bar = ei;
            display('reorthogonal')
        end
    end
    tho = sqrt(x_bar'*fncB(x_bar));% tho = norm(x_bar);
    p = x_bar/tho;
end