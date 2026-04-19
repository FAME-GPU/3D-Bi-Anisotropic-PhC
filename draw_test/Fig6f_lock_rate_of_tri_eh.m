clear; clc; close all
%% ========================= 参数设置 =========================

radius_array = 0.9:0.001:1.3;
problem_all = 'FCC_tri';   % FCC_sphere   FCC_tri   FCC_square  woodpile   sphere_old
problem_short = 'tri'; %   sphere       tri       square    woodpile   sphere_old
gamma_array = 1.009;
n1 = 60; n2 = 60; n3 = 60;
lattice.constant = 3.614910;
lattice.lattice_vector = 0.5*[0,1,1;1,0,1;1,1,0]*lattice.constant;
A = [lattice.lattice_vector(:,1), lattice.lattice_vector(:,2), lattice.lattice_vector(:,3)];  % 3x3
% 固定随机点数量（用于估计“介质/背景”的积分量）
rand_num = 1e5;
% 固定随机种子：保证所有 radius 都用同一批随机点，结果可复现
rng(1);
% 分数坐标（斜坐标）在[0,1)^3里均匀撒点
coef = rand(rand_num, 3);
% 映射到真实空间坐标（用于几何体：内外判断）
points_cart = coef * A;
% 只使用每个分量前两列（i=1,2），分量 k=1..3 => 总共6条曲线
numModesUse = 2;  % i = 1:2
numComp = 3;      % k = 1:3
dim3 = numComp * numModesUse;  % 6
%% ========================= 预计算：插值所需信息 =========================
[i0, j0, k0, u, v, w] = locate_cell_frac_periodic(coef, n1, n2, n3);
W   = trilinear_weights(u, v, w);                                % rand_num x 8
lin = corner_linear_indices_periodic(i0, j0, k0, n1, n2, n3);    % rand_num x 8
%% ========================= 输出数组 =========================
G = numel(gamma_array);
R = numel(radius_array);
lock_rate1_e = zeros(G, R, dim3);  % 介质/全局
lock_rate2_e = zeros(G, R, dim3);  % 介质/背景
lock_rate3_e = zeros(G, R, dim3);  % 平均介质/平均背景
lock_rate1_h = zeros(G, R, dim3);
lock_rate2_h = zeros(G, R, dim3);
lock_rate3_h = zeros(G, R, dim3);
%% ========================= 主循环：gamma -> 读场 -> 插值 -> radius sweep =========================
for gamma_idx = 1:G
    gamma = gamma_array(gamma_idx);
    fprintf("============= 当前 gamma = %.4f =============\n", gamma);
    % filename = sprintf('../Eigen-result/gamma=%.7f.mat', gamma);
    % output = evalc('S=load(filename)');
    % (216000 x 10) 
    filename = sprintf('../test_data_store/%s_n60_gamma_choose4_=%.7f.mat', problem_all, gamma);
    output = evalc('S=load(filename)');
    if iscell(S.result)
        S.result = S.result{1,51};
    end
    e1 = gather(S.result.Eigenmode.e1);
    e2 = gather(S.result.Eigenmode.e2);
    e3 = gather(S.result.Eigenmode.e3);
    h1 = gather(S.result.Eigenmode.h1);
    h2 = gather(S.result.Eigenmode.h2);
    h3 = gather(S.result.Eigenmode.h3);
    % 逐分量(k=1..3)与模态(i=1..2)计算锁光比
    for k = 1:numComp
        for i = 1:numModesUse
            mode_idx = 2*(k-1) + i;   % 第3维索引
            % 取出当前分量、当前模态在网格上的场值（长度=216000）
            switch k
                case 1
                    ve_grid = e1(:,i);  vh_grid = h1(:,i);
                case 2
                    ve_grid = e2(:,i);  vh_grid = h2(:,i);
                case 3
                    ve_grid = e3(:,i);  vh_grid = h3(:,i);
            end
            ve_grid=abs(ve_grid);   vh_grid=abs(vh_grid);
            % ---------- 插值：网格(216000) -> 随机点(rand_num) ----------
            ve_p = interpolate_on_periodic_grid(ve_grid, lin, W);  % rand_num x 1
            vh_p = interpolate_on_periodic_grid(vh_grid, lin, W);  % rand_num x 1
            % 加上虚量
            abs_ve = abs(ve_p);
            abs_vh = abs(vh_p);
            % 全局总量（与 radius 无关）
            tmp1_e = sum(abs_ve);
            tmp1_h = sum(abs_vh);
            % ---------- radius sweep：同一批随机点下统计介质/背景 ----------
            if 1
                r0 = 1.0;
                c  = [0.5 0.5 0.5];      % 缩放中心(分数坐标)
                
                para1_0 = sqrt(3)/2*r0*0.5; 
                para2_0 = 0.5*r0*0.5;
                para3_0 = r0*0.5;
                
                h1_0 = r0 * 0.7;
                h2_0 = (2 - r0) * 0.3;
                
                tri_fru_top0 = [ 0.5-para1_0  0.5-para2_0  h1_0;
                                 0.5+para1_0  0.5-para2_0  h1_0;
                                 0.5          0.5+para3_0  h1_0;
                                 0.5          0.5          0.5;
                                 0.5          0.5          0.5;
                                 0.5          0.5          0.5 ];
                
                tri_fru_bot0 = [ 0.5          0.5          0.5;
                                 0.5          0.5          0.5;
                                 0.5          0.5          0.5;
                                 0.5-para1_0  0.5-para2_0  h2_0;
                                 0.5+para1_0  0.5-para2_0  h2_0;
                                 0.5          0.5+para3_0  h2_0 ];
                center1 = [0.5, 0.5, (h1_0 + 0.5)/2];
                center2 = [0.5, 0.5, (h2_0 + 0.5)/2];
            end
            for rIdx = 1:R
                radius_rate = radius_array(rIdx);
                
                % 几何信息
                s = radius_rate;
                tri_fru_top = [(tri_fru_top0(1:3,:) - center1) .* s + center1; (tri_fru_top0(4:6,:) - center2) .* s + center2];
                tri_fru_bot = [(tri_fru_bot0(1:3,:) - center1) .* s + center1; (tri_fru_bot0(4:6,:) - center2) .* s + center2];
                geometry.tri_fru_top = tri_fru_top * A;
                geometry.tri_fru_bot = tri_fru_bot * A;
                
                % 判断哪些随机点在介质内
                idx_in = geometry_judge(points_cart, geometry);
                
                n_in  = nnz(idx_in);
                n_out = rand_num - n_in;
                % ---- 电场：介质/全局, 介质/背景, 平均介质/平均背景 ----
                tmp2_e = sum(abs_ve(idx_in));    % 介质
                tmp3_e = tmp1_e - tmp2_e;        % 背景
                lock_rate1_e(gamma_idx, rIdx, mode_idx) = safe_div(tmp2_e, tmp1_e);
                lock_rate2_e(gamma_idx, rIdx, mode_idx) = safe_div(tmp2_e, tmp3_e);
                lock_rate3_e(gamma_idx, rIdx, mode_idx) = safe_div(tmp2_e/max(n_in,1), tmp3_e/max(n_out,1));
                % ---- 磁场：同理 ----
                tmp2_h = sum(abs_vh(idx_in));
                tmp3_h = tmp1_h - tmp2_h;
                lock_rate1_h(gamma_idx, rIdx, mode_idx) = safe_div(tmp2_h, tmp1_h);
                lock_rate2_h(gamma_idx, rIdx, mode_idx) = safe_div(tmp2_h, tmp3_h);
                lock_rate3_h(gamma_idx, rIdx, mode_idx) = safe_div(tmp2_h/max(n_in,1), tmp3_h/max(n_out,1));
            end
        end
    end
end
%% ========================= 绘图 (修改部分) =========================
output_folder = '../图-锁光比-插值版';
if ~exist(output_folder, 'dir'); mkdir(output_folder); end
islegend = 0; % 可选：是否显示gamma图例

% 调用新的绘图函数：同时画 E (上排) 和 H (下排)，且仅画第1个特征值
draw_combined_lockrate_mode1(radius_array, gamma_array, lock_rate1_e, lock_rate1_h, '介质-全局', output_folder, islegend);
draw_combined_lockrate_mode1(radius_array, gamma_array, lock_rate2_e, lock_rate2_h, '介质-背景', output_folder, islegend);
draw_combined_lockrate_mode1(radius_array, gamma_array, lock_rate3_e, lock_rate3_h, '平均介质-平均背景', output_folder, islegend);

disp('计算与绘图完成。');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ========================= 函数区 =========================
function draw_combined_lockrate_mode1(radius_array, gamma_array, lock_rate_e, lock_rate_h, suffix, output_folder, islegend)
% 将电场和磁场的3个分量（共6条曲线）合并到同一张图中
% 强制使用对数 Y 轴
% 红色表示电场 (E)，蓝色表示磁场 (H)
% 实线、虚线、点划线分别表示分量 1, 2, 3
    
    fig = figure('Visible','on');
    set(fig, 'Position', [100, 100, 600, 450]); 
    
    myFont = 'Times New Roman';
    labelFontSize = 20;
    axisFontSize = 14;
    
    hold on; box on;
    set(gca, 'FontSize', axisFontSize, 'FontName', myFont); 
    
    lineStyles = {'--', '-', ':'}; 
    
    h_plots = gobjects(0); 
    legend_str = {};       
    
    for k = 1:3
        idx = 2*(k-1) + 1; 
        
        for g = 1:numel(gamma_array)
            % 这里用普通的 plot 就可以，后面统一改坐标轴
            pe = plot(radius_array, squeeze(lock_rate_e(g, :, idx)), ...
                'Color',  [81, 188, 133]/255, 'LineStyle', lineStyles{k}, 'LineWidth', 2);
            
            ph = plot(radius_array, squeeze(lock_rate_h(g, :, idx)), ...
                'Color', [212, 106, 93]/255, 'LineStyle', lineStyles{k}, 'LineWidth', 2);
            
            if g == 1
                h_plots(end+1) = pe;
                legend_str{end+1} = sprintf('$\\mathbf{e}_%d$', k);
                
                h_plots(end+1) = ph;
                legend_str{end+1} = sprintf('$\\mathbf{h}_%d$', k);
            end
        end
    end
    
    % --- 关键修改：强制将 Y 轴设置为对数刻度 ---
    set(gca, 'YScale', 'log');
    
    xlabel('\rho', 'Interpreter','tex', 'FontSize', labelFontSize, 'FontName', myFont);
    ylabel('light localization ratio', 'FontSize', labelFontSize, 'FontName', myFont);
    % title(suffix, 'Interpreter', 'none', 'FontSize', labelFontSize, 'FontName', myFont);
    grid on; 
    
    legend(h_plots, legend_str, 'Interpreter','latex', 'FontSize', 12, 'Location', 'best');
    
    safe_suffix = regexprep(suffix, '[ /\\\(\)]', '_'); 
    safe_suffix = regexprep(safe_suffix, '_+', '_'); 
    if startsWith(safe_suffix, '_'), safe_suffix = safe_suffix(2:end); end
    if endsWith(safe_suffix, '_'), safe_suffix = safe_suffix(1:end-1); end
    
    fname = sprintf('LockingRatio_Mode1_Combined_LogY_%s.png', safe_suffix);
    output_filename = fullfile(output_folder, fname);
    
    exportgraphics(fig, output_filename, 'Resolution', 300);
    % close(fig);
end

function mask = geometry_judge(point, geometry)
    N = size(point,1);
    mask = false(N,1);
    if isfield(geometry, 'quad_fru_top')
        mask = mask | quad_fru_vec(point, geometry.quad_fru_top, geometry.quad_fru_bot);
    end
    if isfield(geometry, 'tri_fru_top')
        mask = mask | tri_fru_vec(point, geometry.tri_fru_top, geometry.tri_fru_bot);
    end
    if isfield(geometry, 'sphere_center')
        mask = mask | sph_judge_vec(point, geometry.sphere_center, geometry.sphere_radius);
    end
end
function inside = quad_fru_vec(points, top, bot, tol)
if nargin < 4, tol = 1e-10; end
N = size(points,1); M = size(top,1) / 4; inside = false(N,1);
for j = 1:M
    T1 = top(4*(j-1)+1,:); T2 = top(4*(j-1)+2,:); T3 = top(4*(j-1)+3,:); T4 = top(4*(j-1)+4,:);
    B1 = bot(4*(j-1)+1,:); B2 = bot(4*(j-1)+2,:); B3 = bot(4*(j-1)+3,:); B4 = bot(4*(j-1)+4,:);
    C = (T1+T2+T3+T4+B1+B2+B3+B4) / 8;
    planesA = [ B1;  T1;  T1;  T2;  T3;  T4 ];
    planesB = [ B2;  T2;  T2;  T3;  T4;  T1 ];
    planesC = [ B3;  T3;  B2;  B3;  B4;  B1 ];
    n = cross(planesB-planesA, planesC-planesA, 2);
    sC = sum((C - planesA) .* n, 2); flip = (sC > 0); n(flip,:) = -n(flip,:);
    ok = true(N,1);
    for m = 1:6
        s = sum((points - planesA(m,:)) .* n(m,:), 2); ok = ok & (s <= tol);
    end
    inside = inside | ok;
end
end
function inside = tri_fru_vec(points, top, bot, tol)
if nargin < 4, tol = 1e-10; end
N = size(points,1); M = size(top,1) / 3; inside = false(N,1);
for j = 1:M
    T1 = top(3*(j-1)+1,:); T2 = top(3*(j-1)+2,:); T3 = top(3*(j-1)+3,:);
    B1 = bot(3*(j-1)+1,:); B2 = bot(3*(j-1)+2,:); B3 = bot(3*(j-1)+3,:);
    C = (T1+T2+T3+B1+B2+B3) / 6;
    planesA = [ B1;  T1;  T1;  T2;  T3;  B1;  B2;  B3 ];
    planesB = [ B2;  T2;  T2;  T3;  T1;  B2;  B3;  B1 ];
    planesC = [ B3;  T3;  B2;  B3;  B1;  T2;  T3;  T1 ];
    n = cross(planesB-planesA, planesC-planesA, 2);
    sC = sum((C - planesA) .* n, 2); flip = (sC > 0); n(flip,:) = -n(flip,:);
    ok = true(N,1);
    for m = 1:8
        s = sum((points - planesA(m,:)) .* n(m,:), 2); ok = ok & (s <= tol);
    end
    inside = inside | ok;
end
end
function y = safe_div(a, b)
if b == 0, y = NaN; else, y = a / b; end
end
function val_p = interpolate_on_periodic_grid(val_grid, lin, W)
f = val_grid(lin); val_p = sum(W .* f, 2);
end
function [i0, j0, k0, u, v, w] = locate_cell_frac_periodic(frac, n1, n2, n3)
s = mod(frac(:,1), 1); t = mod(frac(:,2), 1); q = mod(frac(:,3), 1);
I = floor(s * n1); J = floor(t * n2); K = floor(q * n3);
u = s * n1 - I; v = t * n2 - J; w = q * n3 - K;
i0 = I + 1; j0 = J + 1; k0 = K + 1;
end
function W = trilinear_weights(u, v, w)
W = zeros(numel(u), 8);
W(:,1) = (1-u).*(1-v).*(1-w); W(:,2) = u.*(1-v).*(1-w);
W(:,3) = u.*v.*(1-w);         W(:,4) = (1-u).*v.*(1-w);
W(:,5) = (1-u).*(1-v).*w;     W(:,6) = u.*(1-v).*w;
W(:,7) = u.*v.*w;             W(:,8) = (1-u).*v.*w;
s = sum(W,2); W = W ./ s;
end
function lin = corner_linear_indices_periodic(i0, j0, k0, n1, n2, n3)
i1 = mod(i0, n1) + 1; j1 = mod(j0, n2) + 1; k1 = mod(k0, n3) + 1;
lin = zeros(numel(i0), 8);
lin(:,1) = sub2ind([n1, n2, n3], i0, j0, k0); lin(:,2) = sub2ind([n1, n2, n3], i1, j0, k0);
lin(:,3) = sub2ind([n1, n2, n3], i1, j1, k0); lin(:,4) = sub2ind([n1, n2, n3], i0, j1, k0);
lin(:,5) = sub2ind([n1, n2, n3], i0, j0, k1); lin(:,6) = sub2ind([n1, n2, n3], i1, j0, k1);
lin(:,7) = sub2ind([n1, n2, n3], i1, j1, k1); lin(:,8) = sub2ind([n1, n2, n3], i0, j1, k1);
end
function idx = sph_judge_vec(points, center, r)
dif = points - center; idx = sum(dif.^2, 2) < r^2;
end