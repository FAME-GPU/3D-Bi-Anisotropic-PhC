%% 绘图脚本：5 Solvers + 3 Gammas
clear; clc; close all;

save_filename = 'gamma_matrices_results.mat';
load(save_filename);
dimensions = 4 .* (N_list).^3;
solver_colors = [
    212, 106, 93;  % 红色 bicgstabl
    75, 77,  162;  % 紫色 gmres
    85, 152, 204;  % 蓝色 bicg
    246, 176, 64;  % 黄色 pcg
    81, 188, 133]/255; % 绿色 minres

gamma_markers = {'x', 'o', '+'}; % gamma=1,2,3
target_gammas = [0.2, 0.5, 0.8];
lineWidth = 1.2;
markerSize = 7;

figure('Color', 'w', 'Position', [100, 100, 900, 600]);
hold on;

h_solvers = gobjects(length(solver_names), 1);
h_gammas = gobjects(length(target_gammas), 1);

for g_idx = 1:length(target_gammas)
    real_g_idx = find(gamma_list == target_gammas(g_idx));
    current_matrix = gamma_matrices{real_g_idx};
    
    for s_idx = 1:length(solver_names)
        solver_time = current_matrix(s_idx, :) / 10;
        
        % 绘图
        p = plot(dimensions, solver_time, ...
             'Color', solver_colors(s_idx, :), ...
             'Marker', gamma_markers{g_idx}, ...
             'MarkerSize', markerSize, ...
             'LineWidth', lineWidth, ...
             'HandleVisibility', 'off'); % 隐藏所有原始线条的图例
         
        if g_idx == 2
            set(p, 'MarkerFaceColor', 'none');
        end
    end
end


for s = 1:length(solver_names)
    h_solvers(s) = plot(NaN, NaN, 'Color', solver_colors(s,:), 'LineWidth', 2, ...
        'Marker', 'none', 'DisplayName', solver_names{s});
end

for g = 1:length(target_gammas)
    m_face = 'k'; if g == 2, m_face = 'none'; end % 示意图中用黑色表示标记样式
    h_gammas(g) = plot(NaN, NaN, 'Color', 'k', 'LineStyle', 'none', ...
        'Marker', gamma_markers{g}, 'MarkerSize', markerSize, ...
        'MarkerFaceColor', m_face, 'DisplayName', sprintf('\\gamma = %.1f', target_gammas(g)));
end
all_times = [];
for g_idx = 1:length(target_gammas)
    real_g_idx = find(gamma_list == target_gammas(g_idx));
    current_matrix = gamma_matrices{real_g_idx};
    all_times = [all_times, current_matrix(:)/10];
end

set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'FontName', 'Times New Roman', 'FontSize', 14);

xlim([min(dimensions), max(dimensions)]); 
ylim([min(min(all_times)), max(max(all_times))]);
xlabel('$4n=4n_1n_2n_3$', 'Interpreter', 'latex', 'FontSize', 26,  'FontName', 'Times New Roman');
ylabel('time (s)', 'Interpreter', 'latex', 'FontSize', 26,  'FontName', 'Times New Roman');
% title('The running time of different solvers', 'FontSize', 14);


lgd = legend([h_solvers; h_gammas], 'Location', 'eastoutside');
set(lgd, 'FontSize', 22);
set(get(lgd, 'Title'));


grid on; 
grid minor; 
box on;