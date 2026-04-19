function phys = Physical_settings_sphere(problem_name, gamma, radius)
    global tau
    phys.problem_name = problem_name;
    %% lattice system settings
    phys.lattice.constant = 1;
    if problem_name == "FCC_sphere" || problem_name == "woodpile" || problem_name == "FCC_tri"
        phys.lattice.lattice_vector = 0.5*[0,1,1;1,0,1;1,1,0]*phys.lattice.constant;
    end
    phys.lattice.reciprocal_lattice_vector = 2*pi*inv(phys.lattice.lattice_vector');
    %% geometry settings
    phys.geometry.type = 'Isofunction';
    % isofun: 將滿足 isofun{i}(r,coef) > isoval(i) 的離散點判斷為材料內 (即對應下方的 physical.permittivity{2} 與 physical.permeability{2})
    %   mask: 將滿足 mask{i}(r,coef) == 1 的離散點將強制歸類為背景（即使滿足 isofun{i}(r,coef)>isovol), 注意 mask{i} 僅對 isofun{i} 的判定結果有效
    %      r: 離散點在空間中的位置(r 尺寸須為 Nx3, N 為離散點數目)
    %   coef: 離散點的位置 r 用 unit lattice vector A=[a1,a2,a3] 生成時所對應的係數, 須滿足 r = coef*A'
    a = phys.lattice.constant;
    % 用于计算离散点到多个球心的距离并取最小值，从而判断该点是否在球形区域内
    sphere_isofun = @(r,centers) -min(reshape( vecnorm(kron(ones(size(centers,1),1),r) - kron(centers,ones(size(r,1),1)),2,2), size(r,1), size(centers,1)),[],2);
    
    tau = 0;
    switch phys.problem_name
        case "FCC_sphere"
            a1 = phys.lattice.lattice_vector(:,1);
            a2 = phys.lattice.lattice_vector(:,2);
            a3 = phys.lattice.lattice_vector(:,3);
            phys.geometry.supercell_num = [1,1,1]; 
            sphere_radius = 0.2 * phys.lattice.constant;
             % sphere_radius = sphere_radius * radius;
            sphere_centers = [0.5 0.5 0.5];
            sphere_centers = sphere_centers * (phys.lattice.lattice_vector * diag(phys.geometry.supercell_num))';

            phys.geometry.isoval = -sphere_radius;
            phys.geometry.isofun{1} = @(r,coef) sphere_isofun(r,sphere_centers);            
        case 'woodpile'
            phys.geometry.supercell_num = [1,1,1]; 
            %c = a = 4h;d = √2c = √2a;
            phys.woodpile.width = 0.198*phys.lattice.constant;
            phys.woodpile.height = 1/4*phys.lattice.constant;
            phys.woodpile.distance = phys.lattice.constant/sqrt(2);
            phys.geometry.isofun{1} = @(Point_set) woodpile_isofun(Point_set, phys.woodpile.width, phys.lattice.constant);
    end
    
    
    %% material parameter settings
    phys.permittivity{1} = eye(3);
    phys.permeability{1} = eye(3);    
    phys.magnetoelectric_Xi{1} = zeros(3);    
    phys.magnetoelectric_zeta{1} = zeros(3);
    
    phys.permittivity{2} = diag([12 13 14]);
    phys.permeability{2} = eye(3);
    phys.magnetoelectric_Xi{2} = gamma*1i * eye(3);
    phys.magnetoelectric_zeta{2} = -gamma*1i * eye(3);   
   
    b_inv = inv([phys.permittivity{2} phys.magnetoelectric_Xi{2};phys.magnetoelectric_zeta{2} phys.permeability{2}]);
    b_background_inv = inv([phys.permittivity{1} phys.magnetoelectric_Xi{1};phys.magnetoelectric_zeta{1} phys.permeability{1}]);
    phys.permittivity{2} = b_inv(1:3,1:3);         phys.permittivity{1} = b_background_inv(1:3,1:3);
    phys.magnetoelectric_Xi{2} = b_inv(1:3,4:6);   phys.magnetoelectric_Xi{1} = b_background_inv(1:3,4:6);
    phys.magnetoelectric_zeta{2} = b_inv(4:6,1:3); phys.magnetoelectric_zeta{1} = b_background_inv(4:6,1:3);
    phys.permeability{2} = b_inv(4:6,4:6);         phys.permeability{1} = b_background_inv(4:6,4:6);

    %% boundary condition settings
    phys.BC.x.dir = 'x';
    phys.BC.x.str = 'QP'; % PEC
    phys.BC.y.dir = 'y';
    phys.BC.y.str = 'QP'; % PEC
    phys.BC.z.dir = 'z';
    phys.BC.z.str = 'QP'; % Quasi-Periodic QP
    
    %% Construct Brillouin path
    path_string = 'hGphnGHpn';
    if phys.problem_name == "HEX"
        path_string = 'GKMKHAG';
    elseif phys.problem_name == "FCC" || phys.problem_name == "No227_FCC_SiO2"
       	path_string = 'GXWKGLUWLK|UX';
    elseif phys.problem_name == "No198_SC_FeSi" || phys.problem_name == "SC_test"   % SC
        path_string = 'GXMGRX|MR';
    elseif phys.problem_name == "woodpile"
         path_string = 'XULGXWK';
    end
    phys.lattice.path_string = path_string;

    %% 整理中程式
    B = phys.lattice.reciprocal_lattice_vector;
    
    vertex.G  = B*[   0,   0,   0 ]'; % \gamma
    vertex.N  = B*[   0,   0, 1/2 ]'; % 
    vertex.H  = B*[ 1/2,-1/2, 1/2 ]'; % H
    vertex.P  = B*[ 1/4, 1/4, 1/4 ]'; % 
    vertex.n  = B*[ 1/2,   0,-1/2 ]'; % N
    vertex.p  = B*[ 3/4,-1/4,-1/4 ]'; % P
    vertex.h  = B*[ 1/2, 1/2,-1/2 ]'; % H'
    vertex.q  = 0.5 *vertex.p + 0.5 *vertex.n; % \bar{P}
    vertex.r  = 0.75*vertex.H + 0.25*vertex.G; % \bar{H}
    vertex.m  = 0.2*vertex.H;
    vertex.a  =   vertex.n + vertex.m;
    vertex.b  = - vertex.n + vertex.m;
    vertex.z  = B*[   0,   0,   -1/2 ]';
    vertex.Z  = B*[   0,   0,   1/2 ]';
    
    % fcc
      if phys.problem_name == "FCC" || phys.problem_name == "No227_FCC_SiO2" || phys.problem_name == "woodpile"
        vertex.G = B*[ 0  , 0  , 0   ]';
        vertex.K = B*[ 3/8, 3/8, 3/4 ]';
        vertex.L = B*[ 1/2, 1/2, 1/2 ]';
        vertex.U = B*[ 5/8, 1/4, 5/8 ]';
        vertex.W = B*[ 1/2, 1/4, 3/4 ]';
        vertex.X = B*[ 1/2,   0, 1/2 ]';
      end
    
    % hex
    if phys.problem_name == "HEX"
        % vertex.G  = B*[    0,   0,   0 ]';
        vertex.A  = B*[    0,   0, 1/2 ]';
        vertex.H  = B*[  1/3, 1/3, 1/2 ]';
        vertex.K  = B*[  1/3, 1/3,   0 ]';
        vertex.L  = B*[  1/2,   0, 1/2 ]';
        vertex.M  = B*[  1/2,   0,   0 ]';
    elseif phys.problem_name == "No198_SC_FeSi" || phys.problem_name == "SC_test" %  SC
        vertex.G  = B*[    0,   0,   0 ]';
        vertex.M  = B*[  0.5,  0.5, 0.0 ]';
        vertex.R  = B*[  0.5,  0.5, 0.5 ]';
        vertex.X  = B*[    0,  0.5,   0 ]';
    end

    phys.lattice.vertex = vertex;
    
    n = length(path_string);
    part_num = 10;
    phys.lattice.part_num = part_num;
    wave_vec_array = [];
    for i = 1:n-1
        if strcmp(path_string(i),'|') == 0 && strcmp(path_string(i+1),'|') == 0
            string_x = [ 'subpath(1,:) = linspace( vertex.',path_string(i),'(1), vertex.',path_string(i+1),'(1),', num2str(part_num),');'];
            string_y = [ 'subpath(2,:) = linspace( vertex.',path_string(i),'(2), vertex.',path_string(i+1),'(2),', num2str(part_num),');'];
            string_z = [ 'subpath(3,:) = linspace( vertex.',path_string(i),'(3), vertex.',path_string(i+1),'(3),', num2str(part_num),');'];
            eval(string_x); eval(string_y); eval(string_z);
            if i == n-1 || strcmp(path_string(i+2),'|') == 1
                wave_vec_array = [wave_vec_array,subpath(:,1:end)];
            else
                wave_vec_array = [wave_vec_array,subpath(:,1:end-1)];
            end
        end
    end
    phys.lattice.wave_vec_array = wave_vec_array;
    %% Construct covariant form of material parameters
    A = phys.lattice.lattice_vector;
    phys.permittivity_cov{1} = A'*(phys.permittivity{1}*A);
    phys.permittivity_cov{2} = A'*(phys.permittivity{2}*A);
    phys.permeability_cov{1} = A'*(phys.permeability{1}*A);
    phys.permeability_cov{2} = A'*(phys.permeability{2}*A);    
    
    %% 手性
    phys.magnetoelectric_zeta_cov{1} = A'*(phys.magnetoelectric_zeta{1}*A);
    phys.magnetoelectric_zeta_cov{2} = A'*(phys.magnetoelectric_zeta{2}*A);
    phys.magnetoelectric_Xi_cov{1} = A'*(phys.magnetoelectric_Xi{1}*A);
    phys.magnetoelectric_Xi_cov{2} = A'*(phys.magnetoelectric_Xi{2}*A);
end

function val = multisphere( r,sphere_centers,sphere_radius )
    val = false(size(r,1),1);
    for ii = 1:size(sphere_centers,1)
        val = val | (vecnorm(r-sphere_centers(ii,:),2,2) < sphere_radius) ;
    end
end

function val = cylinder_isofun(r, top, bot)
    num_r = size(r,1);
    num_c = size(top,1);
    kron_r   = kron(ones(num_c,1), r);
    kron_top = kron(top, ones(num_r, 1));
    kron_bot = kron(bot, ones(num_r, 1));
    kron_top_bot = kron(top-bot, ones(num_r,1));
    temp1  =   dot(kron_r - kron_top, kron_top_bot, 2) .* ...
               dot(kron_r - kron_bot, kron_top_bot, 2);
    temp2 =  sqrt( sum(cross(kron_top - kron_r, kron_top_bot, 2).^2, 2)./ sum(kron_top_bot.^2, 2) );
    temp2(temp1>=0) = 1E4;
    val = reshape(temp2, num_r, num_c);
    val = -min(val, [], 2);
end