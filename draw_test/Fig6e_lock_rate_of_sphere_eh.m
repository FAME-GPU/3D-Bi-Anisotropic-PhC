clear; clc; close all
%% ========================= 参数设置 =========================
radius_array = 0.9:0.001:1.3;
problem_all = 'FCC_sphere';   % FCC_sphere   FCC_tri   FCC_square  woodpile   sphere_old
problem_short = 'sphere'; %   sphere       tri       square    woodpile   sphere_old
Gamma_array = 3.46427;
for jj = 1:length(Gamma_array)
gamma_array = Gamma_array(jj);
n1 = 60; n2 = 60; n3 = 60;
lattice.constant = 1;
lattice.lattice_vector = 0.5*[0,1,1;1,0,1;1,1,0]*lattice.constant;
A = [lattice.lattice_vector(:,1), lattice.lattice_vector(:,2), lattice.lattice_vector(:,3)];  % 3x3
rand_num = 1e5;
rng(1);
coef = rand(rand_num, 3);
points_cart = coef * A;
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
    filename = sprintf('../test_data_store/%s_n60_gamma_choose4_=%.7f.mat', problem_all, gamma);
    % filename = sprintf('../Eigen-result/gamma=%.7f.mat', gamma);
    output = evalc('S=load(filename)');
    % (216000 x 10) 
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
            ve_grid=abs(ve_grid);vh_grid=abs(vh_grid);
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
            for rIdx = 1:R
                radius_rate = radius_array(rIdx);
                % 几何信息  正方体
                % % fru_a1 = 0.3*(2-radius_rate); fru_b1 = 0.7*radius_rate;
                % % geometry.quad_fru_top = [ fru_a1  fru_a1  fru_b1;
                % %                                fru_b1  fru_a1  fru_b1;
                % %                                fru_b1  fru_b1  fru_b1;
                % %                                fru_a1  fru_b1  fru_b1];
                % % geometry.quad_fru_bot = [ fru_a1  fru_a1  fru_a1;
                % %                                fru_b1  fru_a1  fru_a1;
                % %                                fru_b1  fru_b1  fru_a1;
                % %                                fru_a1  fru_b1  fru_a1];
                
                % 几何信息  三棱锥
                % % para1 = sqrt(3)/2*radius_rate*0.5; para2 = 0.5*radius_rate*0.5;
                % % para3 = radius_rate*0.5;
                % % h1_point = radius_rate * 0.7; h2_point = (2-radius_rate)*0.3;
                % % geometry.tri_fru_top = [   0.5-para1  0.5-para2   h1_point;
                % %                            0.5+para1  0.5-para2   h1_point;
                % %                            0.5        0.5+para3   h1_point;
                % %                            0.5        0.5         0.5;
                % %                            0.5        0.5         0.5;
                % %                            0.5        0.5         0.5 ] * A;
                % % geometry.tri_fru_bot = [   0.5        0.5         0.5;
                % %                            0.5        0.5         0.5;
                % %                            0.5        0.5         0.5;
                % %                            0.5-para1  0.5-para2   h2_point;
                % %                            0.5+para1  0.5-para2   h2_point;
                % %                            0.5        0.5+para3   h2_point] * A;
                sphere_center = [0.5 0.5 0.5];
                geometry.sphere_center = sphere_center * (lattice.lattice_vector)';
                geometry.sphere_radius = 0.2 * lattice.constant * radius_rate;
                % 判断哪些随机点在介质内
                idx_in = geometry_judge(points_cart, geometry);
                
                % idx_in = sph_judge_vec(points_cart, center_cart, sphere_radius);
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
%% ========================= 绘图 (Modified) =========================
root_output_folder = '../图-锁光比-插值版';
if ~exist(root_output_folder, 'dir'); mkdir(root_output_folder); end
subfolder_name = sprintf('%s-%.7f', problem_all, gamma);
sub_output_folder = fullfile(root_output_folder, subfolder_name);
if ~exist(sub_output_folder, 'dir')
    mkdir(sub_output_folder);
end
islegend = 0; 

draw_combined_lockrate_mode1(radius_array, gamma_array, lock_rate1_e, lock_rate1_h, 'Medium / Total', sub_output_folder, islegend);

draw_combined_lockrate_mode1(radius_array, gamma_array, lock_rate2_e, lock_rate2_h, 'Medium / Background', sub_output_folder, islegend);

draw_combined_lockrate_mode1(radius_array, gamma_array, lock_rate3_e, lock_rate3_h, 'Avg(Med) / Avg(Back)', sub_output_folder, islegend);

disp('计算完成。');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ========================= 函数区 =========================
function draw_combined_lockrate_mode1(radius_array, gamma_array, lock_rate_e, lock_rate_h, suffix, output_folder, islegend)
    
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
    
    % ---  强制将 Y 轴设置为对数刻度 ---
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
        % mask = mask | logical( quad_fru(point, geometry.quad_fru_top, geometry.quad_fru_bot) );
        mask = mask | quad_fru_vec(point, geometry.quad_fru_top, geometry.quad_fru_bot);
    end
    if isfield(geometry, 'tri_fru_top')
        % mask = mask | logical( quad_fru(point, geometry.quad_fru_top, geometry.quad_fru_bot) );
        mask = mask | tri_fru_vec(point, geometry.tri_fru_top, geometry.tri_fru_bot);
    end
    if isfield(geometry, 'sphere_center')
        mask = mask | sph_judge_vec(point, geometry.sphere_center, geometry.sphere_radius);
    end
end
function idx = quad_fru(coef, top, bot)
    n = size(top, 1) / 4;  % n sets of array
    idx = zeros(size(coef,1),1);
    for i = 1:size(coef,1)
        point = coef(i,:);
        for j = 1 : n
            top_1 = top(4*(j-1) + 1, :); bot_1 = bot(4*(j-1) + 1, :);
            top_2 = top(4*(j-1) + 2, :); bot_2 = bot(4*(j-1) + 2, :);
            top_3 = top(4*(j-1) + 3, :); bot_3 = bot(4*(j-1) + 3, :);
            top_4 = top(4*(j-1) + 4, :); bot_4 = bot(4*(j-1) + 4, :);
            % 进行判断
            if (pointSideOfPlane(point, bot_1, top_1, top_2, top_3))
                if (pointSideOfPlane(point, top_1, bot_1, bot_2, bot_3))
                    if (pointSideOfPlane(point, top_4, top_1, top_2, bot_1))
                        if (pointSideOfPlane(point, top_1, top_2, top_3, bot_2))
                            if (pointSideOfPlane(point, top_2, top_3, top_4, bot_3))
                                if (pointSideOfPlane(point, top_3, top_4, top_1, bot_4))
                                    idx(i) = 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
function inside = quad_fru_vec(points, top, bot, tol)
% points: N×3（笛卡尔坐标）
% top/bot: (4*M)×3，M个四棱台，每个四棱台4个顶点（上下底面顺序一致）
% tol: 容差（可选）
if nargin < 4
    tol = 1e-10;
end
N = size(points,1);
M = size(top,1) / 4;
inside = false(N,1);
for j = 1:M
    % 取第 j 个四棱台的 8 个点
    T1 = top(4*(j-1)+1,:);  T2 = top(4*(j-1)+2,:);  T3 = top(4*(j-1)+3,:);  T4 = top(4*(j-1)+4,:);
    B1 = bot(4*(j-1)+1,:);  B2 = bot(4*(j-1)+2,:);  B3 = bot(4*(j-1)+3,:);  B4 = bot(4*(j-1)+4,:);
    % 四棱台中心点（保证在内部，用于校正法向方向）
    C = (T1+T2+T3+T4+B1+B2+B3+B4) / 8;
    % 6 个平面：每个平面用 (a,b,c) 三点定义
    % 注意：这里只要平面覆盖体边界即可；顺序不重要，因为后面会用中心点自动翻转法向
    planesA = [ B1;  T1;  T1;  T2;  T3;  T4 ];  % 每行一个平面的参考点 a
    planesB = [ B2;  T2;  T2;  T3;  T4;  T1 ];  % b
    planesC = [ B3;  T3;  B2;  B3;  B4;  B1 ];  % c
    % 对应含义：
    % 1: bottom 面 (B1,B2,B3)
    % 2: top    面 (T1,T2,T3)
    % 3-6: 四个侧面，大致对应 (T1,T2,B2), (T2,T3,B3), (T3,T4,B4), (T4,T1,B1)
    % 只要是各自的侧面平面即可，中心点校正会保证方向一致。
    % 计算每个平面的法向 n = (b-a)×(c-a)
    n = cross(planesB-planesA, planesC-planesA, 2);   % 6×3
    % 用中心点 C 统一法向方向：要求中心点在“内侧”
    % 规定内侧条件为 (x-a)·n <= 0
    sC = sum((C - planesA) .* n, 2);                  % 6×1
    flip = (sC > 0);                                  % 如果中心点在正侧，则翻转法向
    n(flip,:) = -n(flip,:);
    % 对所有点同时计算 6 个面的半空间判定：
    % s = (p - a)·n ，要求 s <= tol
    % 这里用矩阵化写法：每个面单独算一次（6次），每次是 N×3 点乘
    ok = true(N,1);
    for m = 1:6
        s = sum((points - planesA(m,:)) .* n(m,:), 2); % N×1
        ok = ok & (s <= tol);
        % 小优化：如果已经全 false，可提前 break（可选）
        % if ~any(ok), break; end
    end
    inside = inside | ok; % 多个四棱台取并集
end
end
function inside = tri_fru_vec(points, top, bot, tol)
if nargin < 4
    tol = 1e-10;
end
N = size(points,1);
M = size(top,1) / 3;
inside = false(N,1);
for j = 1:M
    % 取第 j 个三棱台的 6 个点
    T1 = top(3*(j-1)+1,:);  T2 = top(3*(j-1)+2,:);  T3 = top(3*(j-1)+3,:);
    B1 = bot(3*(j-1)+1,:);  B2 = bot(3*(j-1)+2,:);  B3 = bot(3*(j-1)+3,:);
    % 体心点（保证在内部，用于校正法向方向）
    C = (T1+T2+T3+B1+B2+B3) / 6;
    % 5 个面：bottom、top、三块侧面
    % 每个面用 (a,b,c) 三点定义平面
    planesA = [ B1;  T1;  T1;  T2;  T3;  B1;  B2;  B3 ];
    planesB = [ B2;  T2;  T2;  T3;  T1;  B2;  B3;  B1 ];
    planesC = [ B3;  T3;  B2;  B3;  B1;  T2;  T3;  T1 ];
    % 含义：
    % 1: bottom 面 (B1,B2,B3)
    % 2: top    面 (T1,T2,T3)
    % 3: 侧面1  (T1,T2,B2) ~ 连接边 T1-T2 与 B1-B2
    % 4: 侧面2  (T2,T3,B3) ~ 连接边 T2-T3 与 B2-B3
    % 5: 侧面3  (T3,T1,B1) ~ 连接边 T3-T1 与 B3-B1
    % 法向：n = (b-a)×(c-a)
    n = cross(planesB-planesA, planesC-planesA, 2);   % 5×3
    % 统一法向方向：要求中心点在“内侧”
    % 规定内侧条件为 (x-a)·n <= 0
    sC = sum((C - planesA) .* n, 2);                  % 5×1
    flip = (sC > 0);
    n(flip,:) = -n(flip,:);
    % 对所有点一次性判 5 个半空间 (现改为8个，因为考虑到三棱锥的特殊情形)
    ok = true(N,1);
    for m = 1:8
        s = sum((points - planesA(m,:)) .* n(m,:), 2); % N×1
        ok = ok & (s <= tol);
    end
    inside = inside | ok; % 多个三棱台取并集
end
end
function inside = tri_fru_vec_up(points, top, bot, tol) % 上三角的三棱锥
% tri_fru_vec: 判断点是否在三棱台/三棱锥内部（向量化）
%
% points: N×3（笛卡尔坐标）
% top/bot: (3*M)×3，M个三棱台，每个三棱台 top 3点 + bot 3点
% tol: 容差
%
% 返回：
% inside: N×1 logical，true表示在任意一个三棱台/锥内
if nargin < 4
    tol = 1e-10;
end
N = size(points,1);
M = size(top,1) / 3;
inside = false(N,1);
for j = 1:M
    % 取第 j 个三棱台的 6 个点
    T1 = top(3*(j-1)+1,:);  T2 = top(3*(j-1)+2,:);  T3 = top(3*(j-1)+3,:);
    B1 = bot(3*(j-1)+1,:);  B2 = bot(3*(j-1)+2,:);  B3 = bot(3*(j-1)+3,:);
    % 体心点（保证在内部，用于校正法向方向）
    C = (T1+T2+T3+B1+B2+B3) / 6;
    % 5 个面：bottom、top、三块侧面
    % 每个面用 (a,b,c) 三点定义平面
    planesA = [ B1;  T1;  B1;  B2;  B3 ];
    planesB = [ B2;  T2;  B2;  B3;  B1 ];
    planesC = [ B3;  T3;  T2;  T3;  T1 ];
    % 含义：
    % 1: bottom 面 (B1,B2,B3)
    % 2: top    面 (T1,T2,T3)
    % 3: 侧面1  (T1,T2,B2) ~ 连接边 T1-T2 与 B1-B2
    % 4: 侧面2  (T2,T3,B3) ~ 连接边 T2-T3 与 B2-B3
    % 5: 侧面3  (T3,T1,B1) ~ 连接边 T3-T1 与 B3-B1
    % 法向：n = (b-a)×(c-a)
    n = cross(planesB-planesA, planesC-planesA, 2);   % 5×3
    % 统一法向方向：要求中心点在“内侧”
    % 规定内侧条件为 (x-a)·n <= 0
    sC = sum((C - planesA) .* n, 2);                  % 5×1
    flip = (sC > 0);
    n(flip,:) = -n(flip,:);
    % 对所有点一次性判 5 个半空间
    ok = true(N,1);
    for m = 1:5
        s = sum((points - planesA(m,:)) .* n(m,:), 2); % N×1
        ok = ok & (s <= tol);
    end
    inside = inside | ok; % 多个三棱台取并集
end
end
function inside = tri_fru_vec_down(points, top, bot, tol) % 下三角的三棱锥
% tri_fru_vec: 判断点是否在三棱台/三棱锥内部（向量化）
%
% points: N×3（笛卡尔坐标）
% top/bot: (3*M)×3，M个三棱台，每个三棱台 top 3点 + bot 3点
% tol: 容差
%
% 返回：
% inside: N×1 logical，true表示在任意一个三棱台/锥内
if nargin < 4
    tol = 1e-10;
end
N = size(points,1);
M = size(top,1) / 3;
inside = false(N,1);
for j = 1:M
    % 取第 j 个三棱台的 6 个点
    T1 = top(3*(j-1)+1,:);  T2 = top(3*(j-1)+2,:);  T3 = top(3*(j-1)+3,:);
    B1 = bot(3*(j-1)+1,:);  B2 = bot(3*(j-1)+2,:);  B3 = bot(3*(j-1)+3,:);
    % 体心点（保证在内部，用于校正法向方向）
    C = (T1+T2+T3+B1+B2+B3) / 6;
    % 5 个面：bottom、top、三块侧面
    % 每个面用 (a,b,c) 三点定义平面
    planesA = [ B1;  T1;  T1;  T2;  T3 ];
    planesB = [ B2;  T2;  T2;  T3;  T1 ];
    planesC = [ B3;  T3;  B2;  B3;  B1 ];
    % 含义：
    % 1: bottom 面 (B1,B2,B3)
    % 2: top    面 (T1,T2,T3)
    % 3: 侧面1  (T1,T2,B2) ~ 连接边 T1-T2 与 B1-B2
    % 4: 侧面2  (T2,T3,B3) ~ 连接边 T2-T3 与 B2-B3
    % 5: 侧面3  (T3,T1,B1) ~ 连接边 T3-T1 与 B3-B1
    % 法向：n = (b-a)×(c-a)
    n = cross(planesB-planesA, planesC-planesA, 2);   % 5×3
    % 统一法向方向：要求中心点在“内侧”
    % 规定内侧条件为 (x-a)·n <= 0
    sC = sum((C - planesA) .* n, 2);                  % 5×1
    flip = (sC > 0);
    n(flip,:) = -n(flip,:);
    % 对所有点一次性判 5 个半空间
    ok = true(N,1);
    for m = 1:5
        s = sum((points - planesA(m,:)) .* n(m,:), 2); % N×1
        ok = ok & (s <= tol);
    end
    inside = inside | ok; % 多个三棱台取并集
end
end
function bool = pointSideOfPlane(p, q, a, b, c) % point p and q is or not the same side of surface(a,b,c)
    % point = [x,y,z]
    tol = 1e-8;
    pa = a - p;
    qa = a - q;
    normal = cross(b - a, c - a);
    dot_p = dot(pa, normal);
    dot_q = dot(qa, normal);
    if abs(dot_p)< tol || abs(dot_q)< tol
        bool = true; %任一一点在平面上，判定为“同侧”
        return;
    end
    bool = (dot(pa, normal) * dot(qa, normal) > 0);
end
function y = safe_div(a, b)
% 防止除0导致报错
if b == 0
    y = NaN;
else
    y = a / b;
end
end
function val_p = interpolate_on_periodic_grid(val_grid, lin, W)
% 周期网格三线性插值：
% val_grid: (n1*n2*n3) x 1，网格点上的场值
% lin:      rand_num x 8，每个随机点对应8个邻点的线性索引
% W:        rand_num x 8，对应的8个权重
f = val_grid(lin);             % rand_num x 8
val_p = sum(W .* f, 2);        % rand_num x 1
end
function [i0, j0, k0, u, v, w] = locate_cell_frac_periodic(frac, n1, n2, n3)
% 在分数坐标(0~1)空间中定位：
% 网格是 n1*n2*n3 个采样点（周期），相当于把[0,1)均匀采样成 n1 个点
% 对一个点 frac=(s,t,q)，它落在哪个网格立方体由 floor(s*n1) 决定。
s = mod(frac(:,1), 1);
t = mod(frac(:,2), 1);
q = mod(frac(:,3), 1);
I = floor(s * n1);  % 0..n1-1
J = floor(t * n2);
K = floor(q * n3);
% 局部坐标（0~1）
u = s * n1 - I;
v = t * n2 - J;
w = q * n3 - K;
% 转成Matlab 1-based索引（1..n）
i0 = I + 1;
j0 = J + 1;
k0 = K + 1;
end
function W = trilinear_weights(u, v, w)
% 三线性插值权重，对应8个角点顺序：
% [000 100 110 010 001 101 111 011]
W = zeros(numel(u), 8);
W(:,1) = (1-u).*(1-v).*(1-w);
W(:,2) = u.*(1-v).*(1-w);
W(:,3) = u.*v.*(1-w);
W(:,4) = (1-u).*v.*(1-w);
W(:,5) = (1-u).*(1-v).*w;
W(:,6) = u.*(1-v).*w;
W(:,7) = u.*v.*w;
W(:,8) = (1-u).*v.*w;
% 数值归一化，保证和为1
s = sum(W,2);
W = W ./ s;
end
function lin = corner_linear_indices_periodic(i0, j0, k0, n1, n2, n3)
% 对每个点的“左下角网格点”(i0,j0,k0)，找它的8个邻点索引。
% 因为是周期网格，所以 i=n1 的右邻点是 1。
i1 = mod(i0, n1) + 1;
j1 = mod(j0, n2) + 1;
k1 = mod(k0, n3) + 1;
lin = zeros(numel(i0), 8);
lin(:,1) = sub2ind([n1, n2, n3], i0, j0, k0); % 000
lin(:,2) = sub2ind([n1, n2, n3], i1, j0, k0); % 100
lin(:,3) = sub2ind([n1, n2, n3], i1, j1, k0); % 110
lin(:,4) = sub2ind([n1, n2, n3], i0, j1, k0); % 010
lin(:,5) = sub2ind([n1, n2, n3], i0, j0, k1); % 001
lin(:,6) = sub2ind([n1, n2, n3], i1, j0, k1); % 101
lin(:,7) = sub2ind([n1, n2, n3], i1, j1, k1); % 111
lin(:,8) = sub2ind([n1, n2, n3], i0, j1, k1); % 011
end
function idx = sph_judge_vec(points, center, r)
% 判断点是否在球内（真实空间坐标）
dif = points - center;
idx = sum(dif.^2, 2) < r^2;
end