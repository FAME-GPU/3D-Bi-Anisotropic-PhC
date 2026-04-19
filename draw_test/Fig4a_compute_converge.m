clc; clear; close all;
start = 9;
load('ew_converge_result.mat');

ew_converge = real(ew_converge); 

ew_ref = ew_converge(end, :);          
ew_approx = ew_converge(1:end-1, :);   
rel_err = abs(ew_approx - ew_ref) ./ abs(ew_ref);
y_data = log2(rel_err);
x_data = [1/9, 1/18, 1/36, 1/72];

figure('Color', 'w');
hold on;
markers = {'o-', 's-', 'd-', '^-', 'v-', '>-'}; 
colors = lines(6);
legend_str = cell(1, 6); 

for j = 1:6
    p_coeffs = polyfit(log2(x_data), y_data(:, j)', 1);
    convergence_order = p_coeffs(1); 
    legend_str{j} = sprintf('$\\lambda_{%d} \\ (p_{%d} \\approx %.2f)$', j, j, convergence_order);
    
    plot(x_data, y_data(:, j), markers{j}, ...
        'LineWidth', 1.5, ...
        'MarkerSize', 8, ...
        'Color', colors(j,:), ...
        'MarkerFaceColor', colors(j,:));
end
hold off;

ax = gca;
set(ax, 'FontName', 'Times New Roman', 'FontSize', 14); 

set(ax, 'XScale', 'log');     
set(ax, 'XDir', 'reverse');
sorted_ticks = sort(x_data);
xticks(sorted_ticks);         

xticklabels({'1/72', '1/36', '1/18', '1/9'}); 

xlabel('$1/n_1$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
ylabel('$\log_2 \frac{R_m^{2s_1}}{R_m^{s_1}}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
legend(legend_str, 'Interpreter', 'latex', 'Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 15);

grid on;
set(ax, 'GridAlpha', 0.3);
set(ax, 'LineWidth', 1.1);
box on;