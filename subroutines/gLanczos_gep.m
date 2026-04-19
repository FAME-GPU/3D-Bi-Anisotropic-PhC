function [EV,EW,nit] = gLanczos_gep(fncA,fncB,fncinvB,n,eigen_wanted,option)
    m         = option.m;     % size of expanded subspace
    k         = eigen_wanted;  % size of reduced subspace
    nit       = 0;
    itmax     = option.itmax;
    tolerence = option.tol;

    dstart_vec = gpuArray(option.v0);
    dQ         = zeros(  n, m+1, 'gpuArray');
    dalpha     = zeros(  m,   1, 'gpuArray');
    dbeta      = zeros(  m,   1, 'gpuArray');
    dT         = zeros(  m,   m, 'gpuArray');
    de         = zeros(  m,   1, 'gpuArray');
    dU_total   = zeros(  m,   m, 'gpuArray');
    dU         = zeros(  m,   m, 'gpuArray');
    dR         = zeros(  m,   m, 'gpuArray');
    db         = zeros(  m,   1, 'gpuArray');
    dP         = zeros(  n,   1, 'gpuArray');
    dtmp       = zeros(5*n,   1, 'gpuArray');
    dI         =   eye(  m,      'gpuArray');


    dstart_vec = dstart_vec/(dstart_vec'*fncB(dstart_vec));
    
    dQ(1:n,1) = dstart_vec(1:n,1);
    [dQ, dalpha, dbeta] = Lanczos_step_gpu(fncA, fncB, fncinvB, n, dQ, dalpha, dbeta, 0, k, dtmp(1:n), dtmp(n+1:2*n), dtmp(2*n+1:3*n), dtmp(3*n+1:4*n), dtmp(4*n+1:5*n));
    
    while nit < itmax
        nit = nit + 1;
        % expend Krylov subspace by Lanczos
        [dQ, dalpha, dbeta] = Lanczos_step_gpu(fncA, fncB, fncinvB, n, dQ, dalpha, dbeta, k, m, dtmp(1:n), dtmp(n+1:2*n), dtmp(2*n+1:3*n), dtmp(3*n+1:4*n), dtmp(4*n+1:5*n));
        % restart step
        [dQ, dalpha, dbeta] = Implicit_restart_step_gpu(dQ, dalpha, dbeta, k, m, dT, de, dI, dU_total, dU, dR, db, dP);
        % stoping criterion
        abs(dbeta(k))
        while abs(dbeta(k)) < tolerence*abs(dalpha(k))
            T       = diag(dalpha(1:k)) + diag(dbeta(1:k-1), -1) + diag(dbeta(1:k-1), 1);
            [EV,EW] = eigen_solve(dQ, T);
            return
        end
    end
end

%% Lanczos Step
function [dQ, dalpha, dbeta] = Lanczos_step_gpu(fun_A, fncB, fncinvB, n, dQ, dalpha, dbeta, k, m, dw, dr, dtmp_n_1, dtmp_n_2, dtmp_n_3)
    for j = k+1 : m
        dw(1:n,1) = fun_A(dQ(:,j));
        dalpha(j) = real( dQ(:,j)'*dw );
        if j == 1
            dr(1:n,1) = fncinvB(dw(1:n,1)) - dalpha(j)*dQ(1:n,j);
        else
            dr(1:n,1) = fncinvB(dw(1:n,1)) - dalpha(j)*dQ(1:n,j) - dbeta(j-1)*dQ(1:n,j-1);
        end
        [dQ(1:n,j+1), dbeta(j)] = gsreorthog_gpu(dQ(1:n,1:j), fncB, dr(1:n,1), n, j, dtmp_n_1, dtmp_n_2, dtmp_n_3 );
    end
end
%% Implicit Restart step
function [dQ, dalpha, dbeta] = Implicit_restart_step_gpu(dQ, dalpha, dbeta, k, m, dT, de, dI, dU_total, dU, dR, db, dP)
    n = size(dQ,1);
    % construst matrix T
    dT(1:m,1:m) = diag(dalpha(1:m,1)) + diag(dbeta(1:m-1),-1) + diag(dbeta(1:m-1),1);
    % choose unwanted ritz values
    unwanted_ritz = choose_unwanted_gpu(dT,k,m);

    de(1:m-1) = 0; 
    de(m)     = 1;
    
    dU_total(1:m,1:m) = dI(1:m,1:m);
    % single shift step
    for i = 1:m-k
        [dU(1:m,1:m),dR(1:m,1:m)] = qr(dT(1:m,1:m) - unwanted_ritz(i)*dI(1:m,1:m));
        dU_total(1:m,1:m)         = dU_total(1:m,1:m)*dU(1:m,1:m);
        dT(1:m,1:m)               = dR(1:m,1:m)*dU(1:m,1:m) + unwanted_ritz(i)*dI(1:m,1:m);
        for j = 1:m-2
            dT(j,j+2:m) = 0;
        end
        dT(1:m,1:m) = (dT(1:m,1:m) + dT(1:m,1:m)')/2;
    end
    db(1:m)     = dU_total'*de(1:m,1);
    dQ(1:n,1:m) = dQ(1:n,1:m)*dU_total(1:m,1:m);
    
    % new Arnoldi decomposition
    dP(1:n,1)   = dT(k+1,k)*dQ(1:n,k+1) + dbeta(m)*db(k)*dQ(1:n,m+1);
    rho         = norm(dP);
    dQ(1:n,k+1) = dP/rho;
    dbeta(1:k)  = [diag(dT(1:k,1:k),-1);rho];
    dalpha(1:k) = diag(dT(1:k,1:k));
end

function unwanted_ritz = choose_unwanted_gpu(dT,k,m)
    ritz_val      = eig(dT);
    [~,idx]       = sort(abs(ritz_val),'ascend');
    unwanted_ritz = ritz_val(idx(1:m-k));
end

%% Eigen Solver
function [EV,EW] = eigen_solve(Q,T)
    [EV, EW] = eig(T);
    EV       = Q(:,1:length(T))*EV; 
end
%% Reorthogonalize 
function [dp, tho] = gsreorthog_gpu(dQ, fncB, dx, m ,n, dp, dr, dx_bar)
% [p, tho] = gsreorthog(Q, x)
%      x = Q*r + tho*p 
% Orthogonalize x against Q with reorthogonalization
%     r = zeros(n,1);
    tho = dx(1:m,1)'*fncB(dx(1:m,1));
    if tho == 0
        dp(1:n,1) = 0;
        return
    end
    mu = tho;
    dx_bar(1:m,1) = dx(1:m,1);
    
    % find the index of the column with minimal 1-norm
    n1 = dQ(1:m,1)'*fncB(dQ(1:m,1)); %norm(dQ(1:m,1));
    k = 1;
    for i=2:n
        n2 = dQ(1:m,i)'*fncB(dQ(1:m,i));%norm(dQ(1:m,i));
        if(n1>n2)
            n1 = n2; 
            k = i;
        end
    end

    while(1) 
        ds(1:n,1) = dQ(1:m,1:n)'*fncB(dx_bar);
        dr(1:n,1) = dr(1:n,1) + ds(1:n,1);
        % x_bar is the residual
        dx_bar(1:m,1) = dx_bar(1:m,1) - dQ(1:m,1:n)*ds(1:n,1);
        tau = dx_bar(1:m,1)'*fncB(dx_bar(1:m,1)); % norm(dx_bar(1:m,1));
        if(tau/tho >=1/2)
            break;
        end
       
        if(tau > 0.1*mu*eps)
            tho = tau;
        else
            tho = 0.1*tho*eps;
            mu = tho;
            dx_bar(1:m,1) = 0;
            dx_bar(  k,1) = 1;
            display('reorthogonal')
        end
    end
    tho = sqrt(dx_bar(1:m,1)'*fncB(dx_bar(1:m,1))); %norm(dx_bar(1:m,1));
    dp(1:m,1) = dx_bar(1:m,1)/tho;
end