clear; clc;

filename = 'FCC_tri_n60_gamma_choose4_1.009.txt';
% filename = 'FCC_tri_n60_gamma_choose4_0.5.txt';
try
    filetext = fileread(filename);
catch
    error('无法找到或读取文件，请确保文件与此脚本在同一目录下。');
end
pattern = 'rsdl\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*=\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)';

tokens = regexp(filetext, pattern, 'tokens');

num_eigenvalues = 5;
eigen_data = cell(num_eigenvalues, 1);

for k = 1:length(tokens)
    i = str2double(tokens{k}{1}); 
    j = str2double(tokens{k}{2}); 
    residual = str2double(tokens{k}{3});
    if i >= 1 && i <= num_eigenvalues
        eigen_data{i}(1, end+1) = j;
        eigen_data{i}(2, end) = residual;
    end
end

save_struct = struct();

for i = 1:num_eigenvalues
    if ~isempty(eigen_data{i})
        [~, sort_idx] = sort(eigen_data{i}(1, :));
        matrix_to_save = eigen_data{i}(:, sort_idx);
        
        var_name = sprintf('matrix_eig_%d', i);
        save_struct.(var_name) = matrix_to_save;
        
        txt_filename = sprintf('Eigenvalue_%d_Residuals.txt', i);
        writematrix(matrix_to_save, txt_filename, 'Delimiter', '\t');
        
        fprintf('成功生成并保存特征值 %d 的矩阵 -> %s (大小: 2 x %d)\n', i, txt_filename, size(matrix_to_save, 2));
    else
        fprintf('未在文件中找到特征值 %d 的数据。\n', i);
    end
end

save('All_Eigenvalues_Residuals_1p009.mat', '-struct', 'save_struct');
disp('所有非空矩阵已同时打包保存至 All_Eigenvalues_Residuals_1p009.mat');

% save('All_Eigenvalues_Residuals_0p5.mat', '-struct', 'save_struct');
% disp('所有非空矩阵已同时打包保存至 All_Eigenvalues_Residuals_0p5.mat');