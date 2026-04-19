clear; clc; close all;

filename = 'All_Eigenvalues_Residuals_1p009.mat';
% filename = 'All_Eigenvalues_Residuals_0p5.mat';
if ~isfile(filename)
    error('找不到文件: %s，请确保文件位于当前 MATLAB 工作目录下。', filename);
end
data = load(filename);

figure('Name', 'Eigenvalue Residuals', 'Position', [100, 100, 800, 600], 'Color', 'w');
hold on;
box on; 
colors = lines(5); 
markers = {'o', 's', '^', 'd', 'v'};
num_eigenvalues = 5;
for i = 1:num_eigenvalues
    var_name = sprintf('matrix_eig_%d', i); 
    
    if isfield(data, var_name)
        matrix_data = data.(var_name);
        
        iterations = matrix_data(1, :);
        residuals = matrix_data(2, :);
        semilogy(iterations, residuals, '-', ...
                 'Color', colors(i,:), ...
                 'Marker', markers{i}, 'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), ...
                 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('$\\omega_%d$', i));
    end
end

hold off;
ax = gca;
ax.YScale = 'log'; 
power_min = -15; 
power_max = 3;   
ax.YTick = 10.^(power_min : power_max); 
set(gca, 'FontSize', 15, 'FontName', 'Times New Roman');
legend('show', 'Location', 'northeast', 'FontSize', 20, 'Interpreter', 'latex');
xlabel('number of iterations', 'FontSize', 20, 'FontName', 'Times New Roman','Interpreter', 'latex');
ylabel('residual', 'FontSize', 20, 'FontName', 'Times New Roman');

disp('绘图完成！');