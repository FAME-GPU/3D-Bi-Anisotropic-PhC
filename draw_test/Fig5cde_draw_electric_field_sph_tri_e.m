clear; clc; close all;

islabel = true; 
%% ============== 1. 文件夹扫描配置 ==============
task_list = {     
    '..',  'test_data_store',    'sphere_geometry';      
};% 'Sphere',  'Cubes',    'c';      
problem = 'FCC_tri'; % FCC_sphere  FCC_tri  
%% ============== 2. 批量处理主循环 ==============
total_processed = 0;
warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
warning('off', 'MATLAB:dispatcher:load:ReplaceFunctionHandleWithEmpty');
for i = 1:size(task_list, 1)
    main_folder = task_list{i, 1};
    sub_folder  = task_list{i, 2};
    geo_type    = task_list{i, 3};
    
    work_path = fullfile(pwd, main_folder, sub_folder);
    if ~isfolder(work_path), continue; end
    % match_pattern = sprintf('%s_n60_gamma_choose4_=3.4642700.mat', problem);
    match_pattern = sprintf('%s_n60_gamma_choose4_=1.0090000.mat', problem);
    % match_pattern = sprintf('%s_n60_gamma_choose4_=0.5000000.mat', problem);
    files = dir(fullfile(work_path, match_pattern)); 
    if isempty(files), continue; end
    
    fprintf('>>> 处理目录: %s | 几何: %s\n', sub_folder, geo_type);
    
    for k = 1:length(files)
        file_name = files(k).name;
        full_path = fullfile(work_path, file_name);
        save_prefix = sprintf('%s_%s_', main_folder, sub_folder);
            draw_composite_figure(full_path, save_prefix, geo_type, islabel);
            fprintf('    [√] %s\n', file_name);
            total_processed = total_processed + 1;
    end
end
fprintf('\n全部完成！共生成 %d 张图片。\n', total_processed);
%% =========================================================================
%%                             核心绘图函数
%% =========================================================================
function draw_composite_figure(filepath, save_prefix, geometry_type, islabel)
    %% ========== [1. 参数配置] ==========
    % 圆柱
    p_cyl.n_theta = 60; p_cyl.n_z = 20;
    p_cyl.z_bot1=0.264771; p_cyl.z_top1=0.617614;
    p_cyl.z_bot2=p_cyl.z_top1; p_cyl.z_top2=0.735229;
    p_cyl.r_big=0.3321; p_cyl.r_small=0.75*0.3321;
    % 立方体
    % % p_cube.size = 0.49; %gk
    % % p_cube.base1 = [0.01, 0.01, 0.01]; 
    % % p_cube.base2 = [0.99-p_cube.size, 0.99-p_cube.size, 0.99-p_cube.size];
    p_cube.size = 0.4; %gk
    p_cube.base1 = [0.3, 0.3, 0.3]; 
    % 球体
    p_sphere.radius = 0.2; p_sphere.center = [0.5, 0.5, 0.5];
    % 三棱锥
    % 三棱锥(由三棱台退化而来), radius_rate 固定为 1
    p_tri.radius_rate = 1;
    radius_rate = p_tri.radius_rate;
    p_tri.para1 = sqrt(3)/2 * radius_rate * 0.5;
    p_tri.para2 = 0.5       * radius_rate * 0.5;
    p_tri.para3 =             radius_rate * 0.5;
    p_tri.h1    = radius_rate * 0.7;
    p_tri.h2    = (2 - radius_rate) * 0.3;
    % 公共顶点(退化面三个点重合 -> 三棱锥顶点)
    p_tri.apex = [0.5, 0.5, 0.5];
    % 上三棱锥底面(3点)
    p_tri.base_top = [ ...
        0.5 - p_tri.para1,  0.5 - p_tri.para2,  p_tri.h1; ...
        0.5 + p_tri.para1,  0.5 - p_tri.para2,  p_tri.h1; ...
        0.5,                0.5 + p_tri.para3,  p_tri.h1  ...
    ];
    % 下三棱锥底面(3点)
    p_tri.base_bot = [ ...
        0.5 - p_tri.para1,  0.5 - p_tri.para2,  p_tri.h2; ...
        0.5 + p_tri.para1,  0.5 - p_tri.para2,  p_tri.h2; ...
        0.5,                0.5 + p_tri.para3,  p_tri.h2  ...
    ];
    
    %% ========== [数据加载] ==========
    raw = load(filepath);
    if iscell(raw.result)
        result = raw.result{51};
    else
        result = raw.result; 
    end
    opt = result.opt;
    [~, fname, ~] = fileparts(filepath);
    tokens = regexp(fname, 'gamma=([\d\.]+)', 'tokens');
    if ~isempty(tokens), gamma_val = str2double(tokens{1}{1}); else, gamma_val=0; end
    
    A = opt.phys.lattice.lattice_vector; 
    dim = opt.comp.Dim.n1;
    
    % gather_func = @(x) gather(x);
    % gather_func = @(x) x; try; gather_func = @(x) gather(x); catch; end 
    tasks = {
        {gather_func(result.Eigenmode.e1), 'E_x'},
        {gather_func(result.Eigenmode.e2), 'E_y'},
        {gather_func(result.Eigenmode.e3), 'E_z'},
        {gather_func(result.Eigenmode.h1), 'H_x'},
        {gather_func(result.Eigenmode.h2), 'H_y'},
        {gather_func(result.Eigenmode.h3), 'H_z'}
    };
    
    %% ========== [计算统一的包围盒范围] ==========
    %  计算晶格 8 个顶点的最大最小坐标，用于锁定所有子图的视野
    [Vx, Vy, Vz] = ndgrid([0 1]);
    V_lat_corners = [Vx(:), Vy(:), Vz(:)];
    V_cart_corners = (A * V_lat_corners')'; 
    global_xlim = [min(V_cart_corners(:,1)), max(V_cart_corners(:,1))];
    global_ylim = [min(V_cart_corners(:,2)), max(V_cart_corners(:,2))];
    global_zlim = [min(V_cart_corners(:,3)), max(V_cart_corners(:,3))];
    margin = 0.02;
    global_xlim = global_xlim + [-margin, margin];
    global_ylim = global_ylim + [-margin, margin];
    global_zlim = global_zlim + [-margin, margin];

%% ========== [绘图初始化] ==========
alpha_base  = 0.04;  
alpha_scale = 1 - alpha_base; 
alpha_pow   = 1.5;   
max_E = 0; max_H = 0;
for i = 1:3, max_E = max(max_E, max(abs(tasks{i}{1}(:)))); end
for i = 4:6, max_H = max(max_H, max(abs(tasks{i}{1}(:)))); end
if max_E == 0, max_E = 1; end; if max_H == 0, max_H = 1; end
if islabel
    f = figure('Color','w', 'Position', [10, 10, 2000, 900], 'Visible', 'on');
else
    f = figure('Color','w', 'Position', [10, 10, 2000, 900], 'Visible', 'off');
end



right_start = 0.02;      
right_width = 0.90;      
right_height = 0.88;     
col_spacing = -0.4;     
row_spacing = -0.02; 

% 计算每个子图的宽度和高度
tile_width = (right_width - 2*col_spacing) / 3;    % 3列
tile_height = (right_height - row_spacing) / 2;    % 2行

%% --- Tile 2-7: 场分布 ---
grid_lin = linspace(0, 1, dim);
[XX, YY, ZZ] = meshgrid(grid_lin, grid_lin, grid_lin);

field_labels = {
    '$\mathbf{e}_1$', '$\mathbf{e}_2$', '$\mathbf{e}_3$', ...
    '$\mathbf{h}_1$', '$\mathbf{h}_2$', '$\mathbf{h}_3$'
};

% 创建右侧6个子图
for row = 1:1
    for col = 1:3
        % 计算子图索引 (1-6)
        idx = (row-1)*3 + col;
        
        % 计算子图位置
        left = right_start + (col-1)*(tile_width + col_spacing);
        bottom = 0.1 + (2-row)*(tile_height + row_spacing);  % 从下往上计算
        
        % 创建子图
        ax = axes('Position', [left, bottom, tile_width, tile_height]);
        hold(ax, 'on');
        
        % 获取数据
        data = tasks{idx}{1};
        is_E_field = (idx <= 3);
        current_max = is_E_field * max_E + (~is_E_field) * max_H;
        
        % 处理数据
        v_abs = abs(data);
        VV = reshape(v_abs, dim, dim, dim);
        mask = ~isnan(VV); 
        
        P_lat = [XX(mask)'; YY(mask)'; ZZ(mask)'];
        P_real = A * P_lat;
        x_p = P_real(1,:)'; y_p = P_real(2,:)'; z_p = P_real(3,:)';
        v_p = VV(mask);
        
        if ~isempty(x_p)
            v_ratio = v_p / current_max; v_ratio(v_ratio > 1) = 1; 
            alpha_data = alpha_base + alpha_scale * v_ratio.^alpha_pow; 
            alpha_data(alpha_data > 1) = 1;
            
            scatter3(ax, x_p, y_p, z_p, 30, v_p, 'filled', ...
                'MarkerEdgeColor','none', 'MarkerFaceAlpha','flat', 'AlphaData',alpha_data);
        end
        
        colormap(ax, 'parula'); 
        caxis(ax, [0, current_max]);
        
        draw_lattice_box(ax, A); 
        
        view(ax, 30, 30); 
        axis(ax, 'equal', 'off', 'tight');

        xlim(ax, global_xlim); ylim(ax, global_ylim); zlim(ax, global_zlim);
       
        title_pos = [mean(global_xlim), mean(global_ylim), global_zlim(2)*1.05];
        if col == 3
            cb = colorbar(ax); 
            cb.FontName = 'Times New Roman';  
            cb.FontSize = 16;                 
            cb.FontWeight = 'bold';
            cb_pos = cb.Position;
            cb.Position = [ ...
                cb_pos(1) - 0.12, ...           
                cb_pos(2) + cb_pos(4)*0.15, ... 
                cb_pos(3) * 0.8, ...            
                cb_pos(4) * 0.7 ...             
            ];
        end
    end
end

out_name = sprintf('%s%s_FullField_Fixed.png', save_prefix, fname);
exportgraphics(f, out_name, 'Resolution', 300);
if ~islabel
    close(f);
end
end
%% =========================================================================
%%                             通用边框绘制
%% =========================================================================
function draw_lattice_box(ax, A)
    % 风格颜色定义 
    color_dash  = [0.8, 0.6, 0.2]; % 浅橙色 (虚线)
    color_solid = [0.9, 0.7, 0.3]; % 橙色 (实线)
    line_width  = 1.5;             % 统一线宽
    
    % ndgrid 生成顺序:
    [Vx, Vy, Vz] = ndgrid([0 1]);
    V_lat = [Vx(:), Vy(:), Vz(:)];
    V_cart = (A * V_lat')'; 
    
    % 定义所有 12 条边的连接关系
    edges = [
        1 2; 1 3; 1 5; % 从原点出发
        2 4; 2 6;      % 从 a1
        3 4; 3 7;      % 从 a2
        5 6; 5 7;      % 从 a3
        4 8; 6 8; 7 8  % 连向最远点
    ];
    % 绿色框对应的顶点是 3 (0,1,0)
    target_vertex_idx = 3; 
    
    %   创建两个空数组，分别收集虚线和实线
    h_dash = gobjects(0); 
    h_solid = gobjects(0); 
    
    for i = 1:size(edges, 1)
        p1_idx = edges(i, 1);
        p2_idx = edges(i, 2);
        
        pt1 = V_cart(p1_idx, :);
        pt2 = V_cart(p2_idx, :);
        
        % 只要这条边连接到了 3 号点，就画虚线
        if p1_idx == target_vertex_idx || p2_idx == target_vertex_idx
            h = plot3(ax, [pt1(1), pt2(1)], [pt1(2), pt2(2)], [pt1(3), pt2(3)], ...
                'Color', color_dash, 'LineWidth', line_width, 'LineStyle', '--');
            h_dash(end+1) = h; % 放入虚线组
        else
            h = plot3(ax, [pt1(1), pt2(1)], [pt1(2), pt2(2)], [pt1(3), pt2(3)], ...
                'Color', color_solid, 'LineWidth', line_width, 'LineStyle', '-');
            h_solid(end+1) = h; % 放入实线组
        end
    end
    

    % 1. 将虚线强行沉到最底层，让蓝色的场能够盖住它
    if ~isempty(h_dash)
        uistack(h_dash, 'bottom');
    end
    % 2. 将实线强行提到最顶层，防止边缘的散点溢出遮挡它
    if ~isempty(h_solid)
        uistack(h_solid, 'top');
    end
end
%% =========================================================================
%%                             几何体内容绘制
%% =========================================================================
function render_geometry_content(ax, type, A, p_cyl, p_cube, p_sphere, p_tri)
    faces = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];
    
    switch type
        case 'Cylinders'
            [X1, Y1, Z1] = get_cyl_mesh(p_cyl.z_bot1, p_cyl.z_top1, p_cyl.r_big, p_cyl.n_theta, p_cyl.n_z, A);
            surf(ax, X1, Y1, Z1, 'FaceColor','r', 'FaceAlpha', 0.6, 'EdgeColor','none');
            [X2, Y2, Z2] = get_cyl_mesh(p_cyl.z_bot2, p_cyl.z_top2, p_cyl.r_small, p_cyl.n_theta, p_cyl.n_z, A);
            surf(ax, X2, Y2, Z2, 'FaceColor','r', 'FaceAlpha', 0.6, 'EdgeColor','none');
            
        case 'Cubes'
            V1 = get_cube_vertices(p_cube.base1, p_cube.size);
            patch(ax, 'Vertices', (A*V1')', 'Faces', faces, 'FaceColor','r', 'FaceAlpha',0.6, 'EdgeColor','none');
        case 'Sphere'
            [sx, sy, sz] = sphere(60);
            c_real = A * p_sphere.center';
            Xr = c_real(1) + p_sphere.radius * sx;
            Yr = c_real(2) + p_sphere.radius * sy;
            Zr = c_real(3) + p_sphere.radius * sz;
            surf(ax, Xr, Yr, Zr, 'FaceColor','r', 'FaceAlpha',0.6, 'EdgeColor','none');
        case 'Tri'
            % 两个三棱锥: 顶点相同, 底面分别在 z=h1 与 z=h2
            draw_tri_pyramid(ax, A, p_tri.base_top, p_tri.apex);
            draw_tri_pyramid(ax, A, p_tri.base_bot, p_tri.apex);
    end
end
function [X_r, Y_r, Z_r] = get_cyl_mesh(z_bot, z_top, r, n_theta, n_z, A)
    theta = linspace(0, 2*pi, n_theta);
    z = linspace(z_bot, z_top, n_z)';
    [TT, ZZ] = meshgrid(theta, z);
    X_l = 0.5 + r .* cos(TT);
    Y_l = 0.5 + r .* sin(TT);
    Z_l = ZZ;
    coords_r = A * [X_l(:)'; Y_l(:)'; Z_l(:)'];
    X_r = reshape(coords_r(1,:), size(X_l));
    Y_r = reshape(coords_r(2,:), size(Y_l));
    Z_r = reshape(coords_r(3,:), size(Z_l));
end
function V = get_cube_vertices(base, s)
    x0 = base(1); y0 = base(2); z0 = base(3);
    V = [x0 y0 z0; x0+s y0 z0; x0 y0+s z0; x0+s y0+s z0;
         x0 y0 z0+s; x0+s y0 z0+s; x0 y0+s z0+s; x0+s y0+s z0+s];
end
function draw_tri_pyramid(ax, A, base_lat_3x3, apex_lat_1x3)
    V_lat = [base_lat_3x3; apex_lat_1x3];          % 4x3
    V_cart = (A * V_lat')';                        % 4x3, 映射到真实坐标
    faces_tri = [ 1 2 3; 1 2 4; 2 3 4; 3 1 4 ];
    patch(ax, 'Vertices', V_cart, 'Faces', faces_tri, ...
        'FaceColor', 'r', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
end