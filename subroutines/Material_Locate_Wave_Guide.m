function [ Point_idx ] = Material_Locate_Wave_Guide( Point_set, lattice_constant, woodpile_width, airwall_thickness, periodic_num, edge_len )
    Point_idx = 1:size(Point_set,1);
    
    Point_idx_airpile_x = Material_Locate_Airpile( Point_set, lattice_constant, woodpile_width, airwall_thickness, periodic_num(2), periodic_num(3), 2, 3 );    
    Point_idx_airpile_y = Material_Locate_Airpile( Point_set, lattice_constant, woodpile_width, airwall_thickness, periodic_num(3), periodic_num(1), 3, 1 );
    Point_idx_airpile_z = Material_Locate_Airpile( Point_set, lattice_constant, woodpile_width, airwall_thickness, periodic_num(1), periodic_num(2), 1, 2 );
    Point_idx_airwall   = find( (Point_set(:,1) < airwall_thickness(1) ) | (Point_set(:,1) > edge_len(1)-airwall_thickness(1) ) | ...
                                (Point_set(:,2) < airwall_thickness(2) ) | (Point_set(:,2) > edge_len(2)-airwall_thickness(2) ) | ...
                                (Point_set(:,3) < airwall_thickness(3) ) | (Point_set(:,3) > edge_len(3)-airwall_thickness(3) ) );
        
    Point_idx = setdiff(Point_idx, Point_idx_airpile_x);
    Point_idx = setdiff(Point_idx, Point_idx_airpile_y);
    Point_idx = setdiff(Point_idx, Point_idx_airpile_z);
    Point_idx = setdiff(Point_idx, Point_idx_airwall);
end

function [ Point_idx_airpile ] = Material_Locate_Airpile( Point_set, lattice_constant, woodpile_width, airwall_thickness, periodic_num_dir1, periodic_num_dir2, dir1, dir2 )
    
    global flag_plot

    idx_dir1 = []; idx_dir2 = [];
    for i = 1:periodic_num_dir1
        idx_dir1 = union( idx_dir1, ...
                          find( (Point_set(:,dir1) >  airwall_thickness(dir1) + lattice_constant*(i-1) + woodpile_width ) & ...
                                (Point_set(:,dir1) <  airwall_thickness(dir1) + lattice_constant* i    - woodpile_width) ) );
        if flag_plot
            figure(10);subplot(1,3,1);cla
            plot3(Point_set(idx_dir1,1),Point_set(idx_dir1,2),Point_set(idx_dir1,3),'r.');
            axis equal
            axis([min(Point_set(:,1)),max(Point_set(:,1)),min(Point_set(:,2)),max(Point_set(:,2)),min(Point_set(:,3)),max(Point_set(:,3))])
        end
    end
    for i = 1:periodic_num_dir2
        idx_dir2 = union( idx_dir2, ...
                          find( (Point_set(:,dir2) >  airwall_thickness(dir2) + lattice_constant*(i-1) + woodpile_width ) & ...
                                (Point_set(:,dir2) <  airwall_thickness(dir2) + lattice_constant* i    - woodpile_width) ) );
        if flag_plot
            figure(10);subplot(1,3,2);cla
            plot3(Point_set(idx_dir2,1),Point_set(idx_dir2,2),Point_set(idx_dir2,3),'b.');
            axis equal
            axis([min(Point_set(:,1)),max(Point_set(:,1)),min(Point_set(:,2)),max(Point_set(:,2)),min(Point_set(:,3)),max(Point_set(:,3))])
        end
    end
    Point_idx_airpile = intersect( idx_dir1, idx_dir2 );
    if flag_plot
        figure(10);subplot(1,3,3);cla
        plot3(Point_set(Point_idx_airpile,1),Point_set(Point_idx_airpile,2),Point_set(Point_idx_airpile,3),'k.');
        axis equal
        axis([min(Point_set(:,1)),max(Point_set(:,1)),min(Point_set(:,2)),max(Point_set(:,2)),min(Point_set(:,3)),max(Point_set(:,3))])
    end
    
end