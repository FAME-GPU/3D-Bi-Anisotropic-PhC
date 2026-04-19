function [ev, ew, iter] = Lanczos(fncA, fncB, fncBinv, n, k, tol, maxiter, p, v0, flag_gpu)

% v0 = ones(n);
% % v0 = rand(n);
if flag_gpu
    v0 = gpuArray(v0);
end
v0 = v0 / sqrt(abs(v0' * fncB(v0)));
% 
% p = dim;

[ew, ev, iter, res] = Lanczos_iteration(fncA, fncB, fncBinv, n, k, tol, maxiter, p, v0, flag_gpu);

end



function [d, U, kk, res] = Lanczos_iteration(fncA, fncB, fncBinv, n, k, tol, maxiter, p, v, flag_gpu)
% kk 总迭代次数
% mm 重启次数
% Applies Krylov Schur algorithm for Hermitian matrices(reference 1)
% Get previously set algorithm options
eps = 2.2204e-16;

% Initialize variables
Q = zeros(n, p);
if flag_gpu
    Q = gpuArray(Q);
end
% Q[:, 0] = v;
k0 = k;  % k0 is original k. (variable k will be adaptively increased)
Alpha = [];
Beta = [];
alpha = 0;
beta = 0;
d = [];
c = [];
justRestarted = false;
sizeV = 1;
kk = 0;
% fprintf(repmat(' ', 1, 80));  % 63

% Begin Main Algorithm
for mm = 1 : maxiter
    % Build Krylov subspace in Q, H:
    
    for jj = sizeV : p
        
        Q(:, jj) = v;
        r = fncA(Q(:, jj));
        alpha = real(r' * Q(:, jj));
        
        if jj == 1
            pr = fncBinv(r);
            r = pr - alpha * Q(:, jj);
        elseif justRestarted
            pr = fncBinv(r);
            r = GramSchmitGEP(Q, pr, fncB, jj);
            justRestarted = false;
        else
            pr = fncBinv(r);
            r = pr - alpha * Q(:, jj) - beta * Q(:, jj - 1);
            r = r - fncB(Q(:, jj - 1))' * r * Q(:, jj - 1);
        end
        
        r = r - fncB(Q(:, jj))' * r * Q(:, jj);
        
        % Full reorthogonalization
        r = GramSchmitGEP(Q, r, fncB, jj);
        
        
        rr = fncB(r);
        beta = sqrt(r' * rr);
        v = r / beta;
        
        if flag_gpu
            alphac = gather(alpha);
            betac = gather(beta);
        else
            alphac = alpha;
            betac = beta;
        end
                
        Alpha = [Alpha; alphac]; %#ok<AGROW>
        Beta = [Beta; betac];%#ok<AGROW>
        
        kk = kk + 1;
        % fprintf(repmat('\b', 1, 110));  % 65
        % fprintf('Lanczos step: %3d/%3d, restart number: %3d, total_steps: %4d.\n', jj, p, mm - 1, kk);
    end
    
    % Build matrix H
    H1 = diag(Alpha) + diag(Beta(1:end-1), 1) + diag(Beta(1:end-1), -1);
    H = blkdiag(diag(d), H1);
    if ~isempty(d)
        H( 1:k , k+1 ) = c';
        H( k+1 , 1: k ) = c;
    end
    Alpha = [];
    Beta = [];
    
    % Compute eigenvalues
    [U, d] = eig(H, 'vector');
    
    % Implicitly calculate residuals
    res = abs(beta * U(end, :));
    
    % Sort eigenvalues and residuals
    [~, ind] = sort(real(d), 'descend');  % [~, ind] = sort(abs(d), 'descend');
    d = d(ind);
    res = res(ind);
    
    % Number of converged eigenpairs:
    isNotConverged = ~(res(1:k0)' < tol*max(eps^(2/3), abs(d(1:k0))));
    nconv = nnz(~isNotConverged);
    
    if nconv >= k0 || mm == maxiter
        % Stop the algorithm now
        break;
    else
        % Adjust k to prevent stagnating (see reference 2)
        k = k0 + min(nconv, floor((p - k0) / 2));
        if k == 1 && p > 3
            k = floor(p / 2);
        end
    end
    
    
    % Find k most desired eigenvalues of d
    ind = ind(1 : k);
    U = U(:, ind);
    
    % Store variables for next iteration
    if flag_gpu
        Q(:, 1 : k) = Q * gpuArray(U);
    else
        Q(:, 1 : k) = Q * U;
    end
    d = d(1 : k);
    c = betac * U(end, :);
    justRestarted = true;
    sizeV = k + 1;
    % fprintf(nofile,'%d eigenvalues have been caculated.\n', nconv+1);
end

U = U(:, ind(1 : k0));
d = d(1 : k0);

if flag_gpu
    U = Q * gpuArray(U);
else
    U = Q * U;
end
    

end


function r = GramSchmitGEP(V, r, fncB, jj)

for k = 1 : jj
    tmp = fncB(V(:, k))' * r;
    r = r - tmp * V(:, k);

end

end

