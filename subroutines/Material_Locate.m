function [ Edgecenter_inner_idx_x, Edgecenter_inner_idx_y, Edgecenter_inner_idx_z, ...
           Facecenter_inner_idx_x, Facecenter_inner_idx_y, Facecenter_inner_idx_z ] = ...
                Material_Locate( Edgecenter_point_set, Facecenter_point_set, edge_len )
            
Sphere_center = 0.5*edge_len;
Sphere_radius = 0.15;

Edgecenter_inner_idx_x = Material_Locate_Sphere( Edgecenter_point_set.x, Sphere_radius, Sphere_center );
Edgecenter_inner_idx_y = Material_Locate_Sphere( Edgecenter_point_set.y, Sphere_radius, Sphere_center );
Edgecenter_inner_idx_z = Material_Locate_Sphere( Edgecenter_point_set.z, Sphere_radius, Sphere_center );
Facecenter_inner_idx_x = Material_Locate_Sphere( Facecenter_point_set.x, Sphere_radius, Sphere_center );
Facecenter_inner_idx_y = Material_Locate_Sphere( Facecenter_point_set.y, Sphere_radius, Sphere_center );
Facecenter_inner_idx_z = Material_Locate_Sphere( Facecenter_point_set.z, Sphere_radius, Sphere_center );
            
            
end


function [ Point_idx ] = Material_Locate_Sphere( Point_set, Sphere_radius, Sphere_center )
Sum = ( Point_set(:,1) - Sphere_center(1) ).^2 + ( Point_set(:,2) - Sphere_center(2) ).^2 + ( Point_set(:,3) - Sphere_center(3) ).^2;
Point_idx =  find(Sum <= Sphere_radius*Sphere_radius);
end