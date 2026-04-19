function Plot_Eigenmode( Vector, edge_len, grid_num, hax, mode, clim, slice_pt )
x =linspace(0, edge_len(1), grid_num(1));
y =linspace(0, edge_len(2), grid_num(2));
z =linspace(0, edge_len(3), grid_num(3));

% clim = max(abs(real(Vector)));

slice_num = 30;

[X, Y, Z] = meshgrid(x, y, z);

Vector = reshape(Vector, grid_num(1), grid_num(2), grid_num(3));
% Vector(abs(Vector)<1e-4) = 0;

W = zeros(grid_num(2), grid_num(1), grid_num(3));
for t = 1 : grid_num(3)
    W(:, :, t) = Vector(:, :, t)';
end
switch mode
    
    case 1
        %for t = linspace(0, edge_len(1), slice_num)
        for t = slice_pt
            obj = slice(hax, X, Y, Z, W, t, [], [], 'nearest');
            hold(hax, 'on')
            obj.FaceColor = 'interp';
            obj.EdgeColor = 'none';
            alpha(hax, 'color');
            alpha(hax, 'scaled');
        end
    case 2
        %for t = linspace(0, edge_len(2), slice_num)
        for t = slice_pt
            obj = slice(hax, X, Y, Z, W, [], t, [], 'nearest');
            hold(hax, 'on')
            obj.FaceColor = 'interp';
            obj.EdgeColor = 'none';
            alpha(hax, 'color');
            alpha(hax, 'scaled');
        end
    case 3
        %for t = linspace(0, edge_len(3), slice_num)
        for t = slice_pt
            obj = slice(hax, X, Y, Z, W, [], [], t, 'nearest');
            hold(hax, 'on')
            obj.FaceColor = 'interp';
            obj.EdgeColor = 'none';
            alpha(hax, 'color');
            alpha(hax, 'scaled');
            alpha(hax, 0.9);
        end
end

caxis([-clim,clim]);
% colormap hot

obj.DiffuseStrength = 0.8;
colorbar(hax,'horiz')

hold(hax, 'off');
axis(hax, 'equal');
axis(hax, 'tight');

set(hax,'xtick',[])
set(hax,'ytick',[])
set(hax,'ztick',[])
end