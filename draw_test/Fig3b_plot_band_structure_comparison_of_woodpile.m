clear global
clc
close all
clear

load('../test_data_store/woodpile_n90bandstucture_Lanczos_gamma_=0.2000000.mat')
ew = gather(Freq_mtx);
ew = ew/(2*pi).*4;
[numBands, numK] = size(ew);
k = 1:numK; 

data_raw = readmatrix('eig_of_comsol.txt', 'NumHeaderLines', 5);
freq_data_Hz = data_raw(:, end); 
bands_Hz = freq_data_Hz'; 
eig_comsol = reshape(bands_Hz(1:550), 10, []); 

scale_factor = 299792458;
comsol_freq = real(eig_comsol) ./ scale_factor;

[n_modes_c, n_k_c] = size(comsol_freq);
k_comsol = linspace(1, numK, n_k_c);

figure('Color', 'w', 'Position', [100, 100, 600, 450]); 
hold on; 

set(gca, 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.5, 'TickDir', 'in');

k_tick = [1 10 19 28 37 46 55];
k_label = {'X','U','L','\Gamma','X','W','K'}; 
yl = [0, max(ew(:)) * 1.05]; 
for i = 1:length(k_tick)
    xline(k_tick(i), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0, 'HandleVisibility', 'off');
end

band_color = [95/255, 172/255, 228/255];
ew = ew(1:8,:);comsol_freq = comsol_freq(1:8,:);
h1 = plot(k, ew', 'Color', band_color, 'LineWidth', 1.5);

comsol_color = 'r'; 
h2 = plot(k_comsol, comsol_freq, '+', ...
    'Color', comsol_color, ...
    'MarkerSize', 6, ...
    'LineWidth', 1.2);

box on;
xlim([1, numK]);       
set(gca, 'XTick', k_tick);
set(gca, 'XTickLabel', k_label, 'TickLabelInterpreter', 'tex'); 
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14); 
xlabel('wave vector $\mathbf{k}$', 'FontSize', 20, 'FontWeight', 'bold', ...
    'FontName', 'Times New Roman', 'Interpreter', 'latex');
ylabel('frequency $(\omega a/2\pi c)$', 'FontSize', 20, 'FontWeight', 'bold', ...
    'FontName', 'Times New Roman', 'Interpreter', 'latex');
% ylim([0,0.2]);
set(gca, 'Layer', 'top'); 

legend([h1(1), h2(1)], {'our method', 'COMSOL'}, ...
    'Location', 'best', 'Box', 'on', 'Interpreter', 'none', 'FontSize', 14);