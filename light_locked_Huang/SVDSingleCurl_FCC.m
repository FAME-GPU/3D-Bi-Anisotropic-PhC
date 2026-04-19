function [ FFT_parameter, FFT_parameter_bar, SVD_curl ] = SVDSingleCurl_FCC( n_1, n_2, ...
                                    n_3, vec_k, period_vec1, period_vec2, period_vec3, delta )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[mtx_D_kx, mtx_D_jx, mtx_D_jell, Lambda, Lambda_x, Lambda_y, Lambda_z] = construct_fft_mtx_1219(n_1, n_2, ...
                                    n_3, vec_k, period_vec1, period_vec2, period_vec3, delta);

FFT_parameter.n1         = n_1;
FFT_parameter.n2         = n_2;
FFT_parameter.n3         = n_3;
FFT_parameter.mtx_D_kx   = mtx_D_kx;
FFT_parameter.mtx_D_jx   = mtx_D_jx;
FFT_parameter.mtx_D_jell = mtx_D_jell;

FFT_parameter_bar.n1         = n_1;
FFT_parameter_bar.n2         = n_2;
FFT_parameter_bar.n3         = n_3;
FFT_parameter_bar.mtx_D_kx   = sparse(1:n_1,1:n_1,conj(FFT_parameter.mtx_D_kx));
FFT_parameter_bar.mtx_D_jx   = conj(FFT_parameter.mtx_D_jx);
FFT_parameter_bar.mtx_D_jell = conj(FFT_parameter.mtx_D_jell);

FFT_parameter.mtx_Fz     = zeros(n_3, n_1*n_2);
for ii = 1:n_1
    FFT_parameter.mtx_Fz(:,(ii-1)*n_2+1:ii*n_2) = mtx_D_jell(:,:,ii);
end
FFT_parameter.mtx_Fy     = kron(ones(1,n_3), mtx_D_jx);

[Lambda_q, Lambda_Q11, Lambda_Q12]  = construct_Lambda_dat(Lambda_x, Lambda_y, Lambda_z);
[Lambda_q2, Lambda_Q21, Lambda_Q22] = construct_Lambda_dat(-conj(Lambda_x), -conj(Lambda_y), Lambda_z);
[Lambda_q3, Lambda_Q31, Lambda_Q32] = construct_Lambda_dat(-conj(Lambda_x), Lambda_y, -conj(Lambda_z));
[Lambda_q4, Lambda_Q41, Lambda_Q42] = construct_Lambda_dat(Lambda_x, -conj(Lambda_y), -conj(Lambda_z));
 
SVD_curl.Lambda_x         = Lambda_x;
SVD_curl.Lambda_y         = Lambda_y;
SVD_curl.Lambda_z         = Lambda_z;
% SVD_curl.Lambda_Q0        = Lambda_Q0;

n                         = n_1 * n_2 * n_3;
%
[i_idx_1, j_idx_1, val_1] = find(Lambda_Q11);
[i_idx_2, j_idx_2, val_2] = find(Lambda_Q12);

[i_idx_21, j_idx_21, val_21] = find(Lambda_Q21);
[i_idx_22, j_idx_22, val_22] = find(Lambda_Q22);

[i_idx_31, j_idx_31, val_31] = find(Lambda_Q31);
[i_idx_32, j_idx_32, val_32] = find(Lambda_Q32);

[i_idx_41, j_idx_41, val_41] = find(Lambda_Q41);
[i_idx_42, j_idx_42, val_42] = find(Lambda_Q42);


i_idx_Q                   = [i_idx_1; i_idx_2; 3*n+i_idx_21; 3*n+i_idx_22; 6*n+i_idx_31; 6*n+i_idx_32; 9*n+i_idx_41; 9*n+i_idx_42];
j_idx_Q                   = [j_idx_1; n+j_idx_2; 2*n+j_idx_21; 3*n+j_idx_22; 4*n+j_idx_31; 5*n+j_idx_32; 6*n+j_idx_41; 7*n+j_idx_42];
val_Q                     = [val_1; val_2; val_21; val_22; val_31; val_32; val_41; val_42];

i_idx_P                     = [i_idx_2; i_idx_1; 3*n+i_idx_22; 3*n+i_idx_21; 6*n+i_idx_32; 6*n+i_idx_31; 9*n+i_idx_42; 9*n+i_idx_41];
j_idx_P                     = [j_idx_2; n+j_idx_1; 2*n+j_idx_22; 3*n+j_idx_21; 4*n+j_idx_32; 5*n+j_idx_31; 6*n+j_idx_42; 7*n+j_idx_41];
val_P                       = [-conj(val_2); conj(val_1); -conj(val_22); conj(val_21); -conj(val_32); conj(val_31); -conj(val_42); conj(val_41)];

SVD_curl.Lambda_Qr        = sparse([i_idx_Q; 12*n+i_idx_P], [j_idx_Q; 8*n+j_idx_P], [val_Q; val_P]);
SVD_curl.Lambda_QrH       = SVD_curl.Lambda_Qr';

SVD_curl.Lambda_Pr        = sparse([i_idx_P; 12*n+i_idx_Q], [j_idx_P; 8*n+j_idx_Q], [val_P; val_Q]);
SVD_curl.Lambda_PrH       = SVD_curl.Lambda_Pr';

SVD_curl.Sigma_r          = [ Lambda_q.^(0.5) ; Lambda_q.^(0.5) ; ...
                              Lambda_q2.^(0.5); Lambda_q2.^(0.5); ...
                              Lambda_q3.^(0.5); Lambda_q3.^(0.5); ...
                              Lambda_q4.^(0.5); Lambda_q4.^(0.5) ]; %[ inv_Lambda_q_s.^(-1); inv_Lambda_q_s.^(-1) ];

end

% =========================== construct_fft_mtx_1219 ======================
function [mtx_D_kx, mtx_D_jx, mtx_D_jell, Lambda, Lambda_x, Lambda_y, Lambda_z] = construct_fft_mtx_1219(n_1, n_2, ...
                    n_3, vec_k, period_vec1, period_vec2, period_vec3, delta)

theta_k  = 2 * pi * 1i * ( vec_k' * period_vec1 ) / n_1;
mtx_D_kx = exp(theta_k*(0:n_1-1)');

mtx_D_jx  = zeros(n_2,n_1);
for jj = 1:n_1
    phi_jx         = 2 * pi * 1i * ( vec_k' * (period_vec2 - 0.5 * period_vec1) - 0.5 * (jj-1) ) / n_2;
    mtx_D_jx(:,jj) = exp( phi_jx * (0:n_2-1)' );
end

mtx_D_jell = zeros(n_3,n_2,n_1);
delta_p    = period_vec3 - ( period_vec1 + period_vec2 ) / 3;
for jj = 1:n_1
    for ll = 1:n_2
        psi_k               = 2 * pi * 1i * ( vec_k' * delta_p - (jj+ll-2)/3 ) / n_3; 
        mtx_D_jell(:,ll,jj) = exp( psi_k * (0:n_3-1)' );
    end
end
      
    lambda_x = zeros(n_1,1);
    lambda_y = zeros(n_1*n_2,1);
    lambda_z = zeros(n_1*n_2*n_3,1);
        
    for j=1:n_1
        theta_j       = 2 * pi * 1i * ( (j-1) + vec_k' * period_vec1 ) / n_1;
        lambda_x(j,1) = ( exp(theta_j) - 1 ) / delta(1);
    end
    Lambda_x = kron(lambda_x,ones(n_2*n_3,1));
    lambda_x = (conj(lambda_x)).*lambda_x;
    
    nrm_k = norm(vec_k);
    for m=1:n_1
        if ( nrm_k == 0 )
            delta_p2 = period_vec2 - 0.5 * period_vec1;
        else
            delta_p2 = period_vec2 - 0.5 * (period_vec1 + (m-1)/norm(vec_k)^2 * vec_k);
        end
        for n=1:n_2
            phi_n                   = 2 * pi * 1i * ( (n-1) + vec_k' * delta_p2 ) / n_2;
            lambda_y((m-1)*n_2+n,1) = (exp(phi_n) - 1) / delta(2);
        end
    end
    Lambda_y = kron(lambda_y,ones(n_3,1));
    lambda_y = (conj(lambda_y)).*lambda_y;
    
    delta_p = period_vec3 - ( period_vec2 + period_vec1 ) / 3;
    k_delta = vec_k' * delta_p;
    for m=1:n_1
        for n=1:n_2
            for l=1:n_3
                psi_l                                 = 2 * pi * 1i * ( (l-1)+ k_delta - (m+n-2)/3) / n_3;
                lambda_z((m-1)*n_2*n_3+(n-1)*n_3+l,1) = (exp(psi_l) - 1) / delta(3);
            end
        end
    end
    Lambda_z = lambda_z;
    lambda_z = conj(lambda_z).*lambda_z;
    
    Lambda = kron(lambda_x,ones(n_2*n_3,1)) + kron(lambda_y,ones(n_3,1)) + lambda_z;
end

% =========================== construct_Lambda_dat  =======================
function [Lambda_q, Lambda_Q1, Lambda_Q2] = construct_Lambda_dat(Lambda_x, Lambda_y, Lambda_z)

Lambda_q       = conj(Lambda_x).*Lambda_x+conj(Lambda_y).*Lambda_y+conj(Lambda_z).*Lambda_z;
Lambda_p       = Lambda_x + Lambda_y + Lambda_z;
Lambda_p       = Lambda_p.*conj(Lambda_p);
        
Lambda_pq_inv_s = 3 * Lambda_q - Lambda_p; 
%
% Compute Lambda_q2_pq := (3 * Lambda_{q}^2 - Lambda_{q} * Lambda_{p})^{1/2}
Lambda_q2_pq    = Lambda_q.*Lambda_pq_inv_s;
Lambda_q2_pq    = Lambda_q2_pq.^(0.5);
%
% Compute Lambda_pq_inv_s := (3 * Lambda_{q} - Lambda_{p})^{1/2}
Lambda_pq_inv_s = Lambda_pq_inv_s.^(0.5);
%
% Compute inv_Lambda_q_s := (Lambda_x' * Lambda_x + Lambda_y' * Lambda_y + Lambda_z' * Lambda_z)^{-1/2}
% inv_Lambda_q_s  = (real(Lambda_q)).^(-0.5);
% %
% % Compute Lambda_q_s := (Lambda_x' * Lambda_x + Lambda_y' * Lambda_y + Lambda_z' * Lambda_z)^{1/2}
% Lambda_q_s      = (real(Lambda_q))^(0.5);

n       = length(Lambda_z);
idx     = zeros(3*n,1);
jdx     = zeros(3*n,1);
val     = zeros(3*n,1);
%Lambda_Q1 = zeros(3*n,1);
%Lambda_Q2 = zeros(3*n,1);

idx(    1:n,  1) = 1:n;
jdx(    1:n,  1) = 1:n;
val(    1:n,  1) = (conj(Lambda_z)-conj(Lambda_y))./Lambda_pq_inv_s;
idx(  n+1:2*n,1) = n+1:2*n;
jdx(  n+1:2*n,1) = 1:n;
val(  n+1:2*n,1) = (conj(Lambda_x)-conj(Lambda_z))./Lambda_pq_inv_s;
idx(2*n+1:3*n,1) = 2*n+1:3*n;
jdx(2*n+1:3*n,1) = 1:n;
val(2*n+1:3*n,1) = (conj(Lambda_y)-conj(Lambda_x))./Lambda_pq_inv_s;
Lambda_Q2        = sparse(idx, jdx, val);

idx(    1:n,  1) = 1:n;
jdx(    1:n,  1) = 1:n;
val(    1:n,  1) = conj(Lambda_z).*(Lambda_z-Lambda_x)-conj(Lambda_y).*(Lambda_x-Lambda_y);
idx(  n+1:2*n,1) = n+1:2*n;
jdx(  n+1:2*n,1) = 1:n;
val(  n+1:2*n,1) = conj(Lambda_x).*(Lambda_x-Lambda_y)-conj(Lambda_z).*(Lambda_y-Lambda_z);
idx(2*n+1:3*n,1) = 2*n+1:3*n;
jdx(2*n+1:3*n,1) = 1:n;
val(2*n+1:3*n,1) = conj(Lambda_y).*(Lambda_y-Lambda_z)-conj(Lambda_x).*(Lambda_z-Lambda_x);
val(    1:n,  1) = val(    1:n,  1)./Lambda_q2_pq;
val(  n+1:2*n,1) = val(  n+1:2*n,1)./Lambda_q2_pq;
val(2*n+1:3*n,1) = val(2*n+1:3*n,1)./Lambda_q2_pq;
Lambda_Q1        = sparse(idx, jdx, val);

end