function [output_vec] = mtx_prod_vec_A_NFGEP(n, input_vec, SVD_curl, mtx_data, ...
    fun_mtx_TH_prod_vec, fun_mtx_T_prod_vec )

%n          = size(input_vec,1) / 4; %n1 * n2 * n3;
output_vec = zeros(32*n,1);
%
% Compute P_{r} * u
% where P_{r} &= (I_3 \otimes T) [ Lambda_P1 Lambda_P2 ]
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

vec1 = mtx_data.xi_d * tmp1 + tmp2;
if (strcmp(mtx_data.Phi_is_diag,'yes'))
    vec1 = vec1./mtx_data.Phi;
else
    vec1 = mtx_data.Perm_amd * (mtx_data.upper_U \ (mtx_data.Low_L \ (mtx_data.Perm_LU * vec1(mtx_data.Perm_amd_vec,:))));
end
tmp1 = tmp1 + mtx_data.zeta_d * vec1;
% vec1 = - vec1;
% tmp1 = - tmp1;
%
% Compute P_{r}^H * tmp1
% where P_{r} = ( I_3 \otimes T) [ Lambda_P1 Lambda_P2 ]
%
for ii = 1:24 
    tmp1((ii-1)*n+1:ii*n,1) = fun_mtx_TH_prod_vec(tmp1((ii-1)*n+1:ii*n,1));
end 
output_vec(1:16*n,1)     = SVD_curl.Lambda_PrH * tmp1;

%
% Compute Q_{r}^H * vec1
% where Q_{r} = ( I_3 \otimes T) [ Lambda_Q1 Lambda_Q2 ]
%
for ii = 1:24   
    vec1((ii-1)*n+1:ii*n,1) = fun_mtx_TH_prod_vec(vec1((ii-1)*n+1:ii*n,1));
end
output_vec(16*n+1:32*n,1)     = SVD_curl.Lambda_QrH * vec1;

