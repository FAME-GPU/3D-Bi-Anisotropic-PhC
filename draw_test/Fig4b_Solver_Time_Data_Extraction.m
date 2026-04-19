%% 提取不同gamma和N对应的求解器时间数据
gamma_list = [0.2, 0.5, 0.8];  
N_list = 12:6:96;          
solver_names = {'minres', 'gmres', 'pcg', 'bicg', 'bicgstabl'}; 
data_folder = '../test_data_store';
gamma_matrices = cell(1, length(gamma_list)); 
for g_idx = 1:length(gamma_list)
    current_gamma = gamma_list(g_idx);
    current_matrix = nan(length(solver_names), length(N_list));
    
    for n_idx = 1:length(N_list)
        current_N = N_list(n_idx);
        filename = sprintf('N%g_gamma%g.mat', current_N, current_gamma);
        mat_filename = fullfile(data_folder, filename);
        try
            load(mat_filename);
            for s_idx = 1:length(solver_names)
                solver_name = solver_names{s_idx};
                current_matrix(s_idx, n_idx) = time_minres.(solver_name);
            end
            
            fprintf('成功加载并提取：%s\n', mat_filename);
        catch ME
            fprintf('警告：无法加载/提取 %s，原因：%s\n', mat_filename, ME.message);
        end
    end
    gamma_matrices{g_idx} = current_matrix;
end

mat_gamma0p2 = gamma_matrices{1};
mat_gamma0p5 = gamma_matrices{2};
mat_gamma0p8 = gamma_matrices{3};

save_filename = 'gamma_matrices_results.mat';
save(save_filename, 'gamma_matrices', 'gamma_list', 'N_list', 'solver_names');
fprintf('\n已成功将gamma_matrices保存到当前目录，文件名称：%s\n', save_filename);

fprintf('\n========== 结果汇总 ==========\n');
fprintf('每个矩阵尺寸：%d行（求解器）× %d列（N值）\n', length(solver_names), length(N_list));
fprintf('gamma=0 矩阵前5列数据：\n');
disp(mat_gamma0p2(:, 1:5));