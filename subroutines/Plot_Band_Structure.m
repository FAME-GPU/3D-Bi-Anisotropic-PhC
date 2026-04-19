function obj = Plot_Band_Structure( path_string, part_num, Freq_array, hax )

    part_len = 0.2; % Set width of partition bar
    
    % Set XData    
    part_bar_idx = [];
    vertex_idx   = 1;
    XData        = 0;
    part_bar_num = 0;
    label_symbol_array{1} = path_string(1);
    if path_string(1) == '|'
        error('Invalid path string format');
    end
    for i = 2:length(path_string)
        if (path_string(i-1) ~= '|') && (path_string(i) ~= '|')
            vertex_idx = [ vertex_idx, (vertex_idx(end)+part_num-1)];
            tmp = linspace( XData(end), (XData(end)+1), part_num);
            XData = [XData, tmp(2:end)];
            label_symbol_array{end+1} = path_string(i);
        elseif (path_string(i-1) == '|') && (path_string(i) ~= '|')
            vertex_idx = [ vertex_idx, (vertex_idx(end)+1) ];
            label_symbol_array{end+1} = path_string(i);
        elseif (path_string(i) == '|')
            part_bar_num = part_bar_num + 1;
            part_bar_idx = [part_bar_idx, [vertex_idx(end);vertex_idx(end)+1]];
            XData = [XData, (XData(end)+part_len)];
        end
    end
    % axes setting
    hold( hax, 'on');
    % Plot Band structure
    for i = 1:size(Freq_array,1)
        YData = Freq_array(i,:);
        obj{i} = plot( hax, XData, YData, '.-', 'Color', [0 0.4470 0.7410], 'MarkerSize', 10 );
    end
    % Plot partition bar
    for i = 1:part_bar_num
        x = [XData(part_bar_idx(1,i))+0.05*part_len, XData(part_bar_idx(2,i))-0.05*part_len, XData(part_bar_idx(2,i))-0.05*part_len, XData(part_bar_idx(1,i))+0.05*part_len];
        y = [0,0,1.1*max(Freq_array(:)),1.1*max(Freq_array(:))];
        fill( hax, x, y, [164,162,166]/255, 'LineStyle', 'none' );
    end
    
    % Set label 
    set(hax,'FontSize'  ,                        12 );
    set(hax,'XTick'     ,         XData(vertex_idx) );
    set(hax,'XTickLabel', upper(label_symbol_array) );  
    if length(XData) ~= 1
        axis( hax,[ min(XData),max(XData),0, 1.1*max(Freq_array(:))] );
    end
    grid( hax, 'on');
end
