function [ ev ] = Compute_eigenvector( GEP_mtx_B, n, input_vec, SVD_curl, mtx_data, fun_mtx_T_prod_vec )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%
% Compute P_{r} * v
% where
%    P_{r} = ( I_3 \otimes T) [ Lambda_P1 Lambda_P2 ] 
%

tmp1 = SVD_curl.Lambda_Pr * input_vec(1:16*n,1); 

for ii = 1:24  
    tmp1((ii-1)*n+1:ii*n,1) = fun_mtx_T_prod_vec(tmp1((ii-1)*n+1:ii*n,1));
end
%
% Compute Q_{r} * v
% where
%    Q_{r} = ( I_3 \otimes T) [ Lambda_Q1 Lambda_Q2 ] 
%
tmp2 = SVD_curl.Lambda_Qr * input_vec(16*n+1:32*n,1); 

for ii = 1:24 
    tmp2((ii-1)*n+1:ii*n,1) = fun_mtx_T_prod_vec(tmp2((ii-1)*n+1:ii*n,1));
end

% ev_tmp = GEP_mtx_B \ [tmp1; tmp2 ];

vec1 = mtx_data.xi_d * tmp1 + tmp2;
if (strcmp(mtx_data.Phi_is_diag,'yes'))
    vec1 = vec1./mtx_data.Phi;
else
    vec1 = mtx_data.Perm_amd * (mtx_data.upper_U \ (mtx_data.Low_L \ (mtx_data.Perm_LU * vec1(mtx_data.Perm_amd_vec,:))));
end
%tmp1 = tmp1 + mtx_data.zeta_d * vec1;

ev   = (-1i) * [ -vec1; tmp1 + mtx_data.zeta_d * vec1 ];

% err = norm(ev - ev_tmp)

% ev = GEP_mtx_B \ [tmp1; tmp2 ];
end

