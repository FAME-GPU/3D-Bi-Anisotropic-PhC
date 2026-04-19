function Fig6cd_plot_max_field(shape_type)
    if ~ismember(shape_type, {'sphere', 'tri'})
        error('输入参数错误：必须是 ''sphere'' 或 ''tri''。');
    end

    % =========================================================
    % 根据输入类型加载对应的基础文件和设置路径模板
    % =========================================================
    if strcmp(shape_type, 'tri')
        base_file = '../test_data_store/FCC_tri_n60_gamma_choose4_=1.0090000.mat';
        file_template = '../radius_data/radius_rate_tri=%.3f.mat';
    else % sphere
        base_file = '../test_data_store/FCC_sphere_n60_gamma_choose4_=3.4642700.mat';
        file_template = '../radius_data/radius_rate_sphere=%.3f.mat';
    end

    fprintf('当前选择形状: %s\n', shape_type);
    fprintf('正在加载基础数据: %s...\n', base_file);
    base_data = load(base_file, 'result');
    
    % 提取并 gather e1, e2, e3 和 h1, h2, h3 的第一列数据
    e1 = gather(base_data.result.Eigenmode.e1);
    e2 = gather(base_data.result.Eigenmode.e2);
    e3 = gather(base_data.result.Eigenmode.e3);
    h1 = gather(base_data.result.Eigenmode.h1);
    h2 = gather(base_data.result.Eigenmode.h2);
    h3 = gather(base_data.result.Eigenmode.h3);
    e1 = e1(:,1); e2 = e2(:,1); e3 = e3(:,1); 
    h1 = h1(:,1); h2 = h2(:,1); h3 = h3(:,1); 
    
    % =========================================================
    % 设置介质半径的循环范围
    % =========================================================
    radii = 0.900:0.010:1.300; 
    num_radii = length(radii);
    max_e_arr = zeros(num_radii, 1);
    max_h_arr = zeros(num_radii, 1);
    
    % =========================================================
    % 循环遍历不同半径的索引文件并计算最大模值
    % =========================================================
    fprintf('开始计算介质外的最大模值...\n');
    for i = 1:num_radii
        r = radii(i);
       
        filename = sprintf(file_template, r); 
        
        if ~exist(filename, 'file')
            warning('未找到文件: %s，跳过该点。', filename);
            max_e_arr(i) = NaN;
            max_h_arr(i) = NaN;
            continue;
        end
        
        idx_data = load(filename, 'Standard_inner_idx');
        inner_idx = idx_data.Standard_inner_idx;
        
        e1_out = e1; e1_out(inner_idx) = 0;
        e2_out = e2; e2_out(inner_idx) = 0;
        e3_out = e3; e3_out(inner_idx) = 0;
        
        max_e_components = [max(abs(e1_out)), max(abs(e2_out)), max(abs(e3_out))];
        max_e_arr(i) = max(max_e_components);
        
        h1_out = h1; h1_out(inner_idx) = 0;
        h2_out = h2; h2_out(inner_idx) = 0;
        h3_out = h3; h3_out(inner_idx) = 0;
        
        max_h_components = [max(abs(h1_out)), max(abs(h2_out)), max(abs(h3_out))];
        max_h_arr(i) = max(max_h_components);
    end
    
    % =========================================================
    % 绘制结果折线图
    % =========================================================
    fig_name = sprintf('Maximum Field Magnitude Outside (%s)', upper(shape_type));
    figure('Name', fig_name, 'Position', [100, 100, 700, 500]);
    redColor   = [212 106 93]/255; 
    blueColor  = [95, 172, 228]/255;
    
    % 绘制 e 场模最大值
    plot(radii, max_e_arr, '-o', 'Color', redColor, 'LineWidth', 1.5, ...
        'MarkerFaceColor', redColor, 'MarkerEdgeColor', redColor);
    hold on;
    % 绘制 h 场模最大值
    plot(radii, max_h_arr, '-s', 'Color', blueColor, 'LineWidth', 1.5, ...
        'MarkerFaceColor', blueColor, 'MarkerEdgeColor', blueColor);
    grid on;
    
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
    xlabel('$\rho$','Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
    ylabel('maximal absolute value', 'FontName', 'Times New Roman', 'FontSize', 20);
    legend({'$m_\mathbf{e}$', '$m_\mathbf{h}$'}, 'Interpreter', 'latex', ...
        'Location', 'best', 'FontSize', 20);
 
    fprintf('计算与绘图完成！\n');
end