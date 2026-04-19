function vec_mtx_T_prod_vec = p_mtx_T_prod_vec_FCC(input_vec, FFT_parameter)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n_1         = FFT_parameter.n1;
n_2         = FFT_parameter.n2;
n_3         = FFT_parameter.n3;

mtx_Q_tilde = reshape(input_vec, n_3, n_1*n_2);
mtx_Q_tilde = FFT_parameter.mtx_Fz.*(ifft(mtx_Q_tilde));
 
mtx_Q_y_tilde = zeros(n_2, n_1*n_3);

for ii = 1:n_1
    mtx_Q_y_tilde(:,(ii-1)*n_3+1:ii*n_3) = mtx_Q_tilde(:,(ii-1)*n_2+1:ii*n_2).';
end

mtx_Q_y_tilde = ifft(mtx_Q_y_tilde);

idx_col = zeros(n_1*n_3,1);
for ii = 1:n_3
    idx_col((ii-1)*n_1+1:ii*n_1,1) = (ii:n_3:(n_1-1)*n_3+ii)';
end

mtx_Q_y_tilde = FFT_parameter.mtx_Fy.*mtx_Q_y_tilde(:,idx_col);

mtx_Q_x_tilde = zeros(n_1, n_2*n_3);
for ii = 1:n_3
    mtx_Q_x_tilde(:,(ii-1)*n_2+1:ii*n_2) = mtx_Q_y_tilde(:,(ii-1)*n_1+1:ii*n_1).';
end

% mtx_Q_x_tilde = ifft(mtx_Q_x_tilde);

vec_mtx_T_prod_vec = sqrt(n_1 * n_2 * n_3) * reshape(sparse(1:n_1,1:n_1,FFT_parameter.mtx_D_kx)*ifft(mtx_Q_x_tilde),n_1*n_2*n_3,1);

end

