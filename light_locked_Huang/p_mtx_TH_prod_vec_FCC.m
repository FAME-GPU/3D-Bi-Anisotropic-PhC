function vec_mtx_TH_prod_vec = p_mtx_TH_prod_vec_FCC(input_vec, FFT_parameter_bar)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n_1           = FFT_parameter_bar.n1;
n_2           = FFT_parameter_bar.n2;
n_3           = FFT_parameter_bar.n3;

% mtx_P_x_tilde = fft(sparse(1:n_1,1:n_1,conj(FFT_parameter.mtx_D_kx)) * reshape(input_vec, n_1, n_3*n_2));
mtx_P_x_tilde = fft(FFT_parameter_bar.mtx_D_kx * reshape(input_vec, n_1, n_3*n_2));

mtx_P_y_tilde = zeros(n_2, n_1*n_3);
for ii = 1:n_1
    mtx_P_y_tilde(:,(ii-1)*n_3+1:ii*n_3) = sparse(1:n_2,1:n_2,FFT_parameter_bar.mtx_D_jx(:,ii)) * reshape(mtx_P_x_tilde(ii,:).', n_2, n_3);
%     mtx_P_y_tilde(:,(ii-1)*n_3+1:ii*n_3) = sparse(1:n_2,1:n_2,conj(FFT_parameter.mtx_D_jx(:,ii))) * reshape(mtx_P_x_tilde(ii,:).', n_2, n_3);
end

mtx_P_y_tilde = fft(mtx_P_y_tilde);

mtx_P_z_tilde = zeros(n_3, n_1*n_2);
for ii = 1:n_1
    mtx_P_z_tilde(:,(ii-1)*n_2+1:ii*n_2) = FFT_parameter_bar.mtx_D_jell(:,1:n_2,ii).*(mtx_P_y_tilde(:,(ii-1)*n_3+1:ii*n_3).');
%     mtx_P_z_tilde(:,(ii-1)*n_2+1:ii*n_2) = conj(FFT_parameter_bar.mtx_D_jell(:,1:n_2,ii)).*(mtx_P_y_tilde(:,(ii-1)*n_3+1:ii*n_3).');
end

vec_mtx_TH_prod_vec = reshape(fft(mtx_P_z_tilde),n_1*n_2*n_3,1) / sqrt(n_1 * n_2 * n_3);
end

