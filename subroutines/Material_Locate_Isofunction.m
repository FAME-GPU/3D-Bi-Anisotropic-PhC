function [ Point_idx ] = Material_Locate_Isofunction( coef, lattice_vec, phys)
    % a1 = lattice_vec(:,1) column vector
    geometry = phys.geometry;
    supercell_lattice_vec = lattice_vec * diag(geometry.supercell_num);
    Point_set = coef*supercell_lattice_vec';
    unitcell_coef = Point_set/lattice_vec';
    for i = 1:length(geometry.isofun)
        if strcmp(phys.problem_name, 'woodpile') 
            f = geometry.isofun{i};
            % Point_idx_tmp{i} = f(Point_set);
            Point_idx_tmp{i} = f(Point_set);
        else
            f = geometry.isofun{i};
            v = geometry.isoval(i);
            Point_idx_tmp{i} = find(f(Point_set,unitcell_coef) > v); % f(Point_set,unitcell_coef)即sphere_isofun(Point_set,sphere_centers)
        end
        if isfield(geometry,'mask')
            mask = geometry.mask{i};
            Point_idx_tmp{i} = setdiff( Point_idx_tmp{i}, find(mask(Point_set,unitcell_coef) == 1) );
        end
    end
    
    Point_idx = [];
    for i = 1:length(geometry.isofun)
        Point_idx = union(Point_idx,Point_idx_tmp{i});
    end
    
    % Quadrilateral Frustum
    if phys.problem_name == "FCC_square"
        idx_qua = find( quad_fru(coef, lattice_vec, geometry, geometry.quad_fru_top, geometry.quad_fru_bot) );
        Point_idx = union(Point_idx, idx_qua);
    elseif phys.problem_name == "FCC_tri"
        idx_tri = find( tri_fru(coef, lattice_vec, geometry, geometry.tri_fru_top, geometry.tri_fru_bot) );
        Point_idx = union(Point_idx, idx_tri);
    end
end

function idx = tri_fru(coef, lattice_vec, geometry, top, bot)
    supercell_lattice_vec = lattice_vec * diag(geometry.supercell_num);
    Point_set = coef * supercell_lattice_vec'; % coef: each point is row vector
    unitcell_coef = Point_set/lattice_vec';
    
    top_coef = top * lattice_vec';
    bot_coef = bot * lattice_vec';

    n = max(size(top, 1) / 3, size(bot, 1) / 3); % the number of geometry needed to compare
    idx = zeros(size(Point_set, 1), 1);
    for i = 1:size(Point_set, 1)
        point = Point_set(i, :);
        for j = 1:n
            top_1 = top_coef(3*(j-1)+1, :); bot_1 = bot_coef(3*(j-1)+1, :);
            top_2 = top_coef(3*(j-1)+2, :); bot_2 = bot_coef(3*(j-1)+2, :);
            top_3 = top_coef(3*(j-1)+3, :); bot_3 = bot_coef(3*(j-1)+3, :);
            % judge
            if sum((top_1-top_2).^2) < 1e-6 % case1
                if pointSideOfPlane(point, top_1, bot_1, bot_2, bot_3)
                    if pointSideOfPlane(point, bot_3, bot_1, bot_2, top_1)
                        if pointSideOfPlane(point, bot_1, bot_2, bot_3, top_1)
                            if pointSideOfPlane(point, bot_2, bot_3, bot_1, top_1)
                                idx(i) = 1;
                            end
                        end
                    end
                end
            elseif sum((bot_1-bot_2).^2) < 1e-6 % case2
                if pointSideOfPlane(point, bot_1, top_1, top_2, top_3)
                    if pointSideOfPlane(point, top_3, top_1, top_2, bot_1)
                        if pointSideOfPlane(point, top_1, top_2, top_3, bot_1)
                            if pointSideOfPlane(point, top_2, top_3, top_1, bot_1)
                                idx(i) = 1;
                            end
                        end
                    end
                end
            else % case3 -- general
                if pointSideOfPlane(point, bot_1, top_1, top_2, top_3)
                    if pointSideOfPlane(point, top_1, bot_1, bot_2, bot_3)
                        if pointSideOfPlane(point, top_3, top_1, top_2, bot_1)
                            if pointSideOfPlane(point, top_1, top_2, top_3, bot_2)
                                if pointSideOfPlane(point, top_2, top_3, top_1, bot_3)
                                    idx(i) = i;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

function idx = quad_fru(coef, lattice_vec, geometry, top, bot)
    supercell_lattice_vec = lattice_vec * diag(geometry.supercell_num);
    Point_set = coef * supercell_lattice_vec'; % coef: each point is row vector
    unitcell_coef = Point_set/lattice_vec';
    
    top_coef = top * lattice_vec';    %gk
    bot_coef = bot * lattice_vec';    %gk
    
    n = size(top, 1) / 4;  % n sets of array
    idx = zeros(size(Point_set,1),1);
    for i = 1:size(Point_set,1)
        point = Point_set(i,:);       %gk 
        for j = 1 : n
            top_1 = top_coef(4*(j-1) + 1, :); bot_1 = bot_coef(4*(j-1) + 1, :);
            top_2 = top_coef(4*(j-1) + 2, :); bot_2 = bot_coef(4*(j-1) + 2, :);
            top_3 = top_coef(4*(j-1) + 3, :); bot_3 = bot_coef(4*(j-1) + 3, :);
            top_4 = top_coef(4*(j-1) + 4, :); bot_4 = bot_coef(4*(j-1) + 4, :);
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

% function bool = pointSideOfPlane(p, q, a, b, c) % point p and q is or not the same side of surface(a,b,c)
%     % point = [x,y,z]
%     pa = a - p;
%     qa = a - q;
%     normal = cross(b - a, c - a);
%     bool = (dot(pa, normal) * dot(qa, normal) > 0);
% end
% gk
function bool = pointSideOfPlane(p, q, a, b, c) % point p and q is or not the same side of surface(a,b,c)
    % point = [x,y,z]
    tol = 1e-8;
    pa = a - p;
    qa = a - q;
    normal = cross(b - a, c - a);
    dot_p = dot(pa, normal);  
    dot_q = dot(qa, normal);  
    if abs(dot_p) < tol || abs(dot_q) < tol
       bool = true;  % 任一一点在平面上，判定为“同侧”
       return;
    end
    bool = (dot(pa, normal) * dot(qa, normal) > 0);
end
 