%% 绘制特征值求解时间和迭代次数随波矢路径变化的图
data_file = '../test_data_store/woodpile_n90bandstucture_Lanczos_gamma_=0.2000000.mat';
if ~exist(data_file, 'file')
    error(['找不到文件: ', data_file, '，请检查文件路径和名称是否正确！']);
end
load(data_file);  

solve_time = time_minres.compute_eigs;  

itno_sum = iter_records;

figure('Position', [100, 100, 800, 600]);  
set(groot, 'DefaultTextInterpreter', 'latex');  
set(groot, 'DefaultAxesFontSize', 14);        

x_pos = [1, 10, 19, 28, 37, 46, 55];  
x_labels = {'X', 'U', 'L', '\Gamma', 'X', 'W', 'K'};  

% ======================  绘制第一张图：求解时间 ======================
subplot(2, 1, 1);  
plot(1:55, solve_time', 'o-', 'Color', [81, 188, 133]/255, 'LineWidth', 1.2, 'MarkerSize', 4);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14); 
xlim([1, 55]);  
ylim([50, 90]);  
xticks(x_pos);   
xticklabels(x_labels);  
xlabel('wave vector($\mathbf{k}$)', 'FontName', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('time(s)', 'FontName', 'Times New Roman', 'FontSize', 24);  
grid on;  

% ====================== 绘制第二张图：迭代次数 ======================
subplot(2, 1, 2);  
plot(1:55, itno_sum, 's-', 'Color', [212, 106, 93]/255, 'LineWidth', 1.2, 'MarkerSize', 4);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14); 
xlim([1, 55]);
ylim([90, 150]); 
xticks(x_pos);
xticklabels(x_labels);
xlabel('wave vector($\mathbf{k}$)', 'FontName', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('iter number', 'FontName', 'Times New Roman', 'FontSize', 24);
grid on;  