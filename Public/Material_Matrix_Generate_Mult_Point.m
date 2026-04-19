function [ mtx_B ] = Material_Matrix_Generate_Mult_Point( opt ,type , radius_rate)
% 3*3分块对角矩阵M,N,L,J,M_inv??phi_inv??
% 注释31 32行 不生成L 2025.10.30 
    Dim = opt.comp.Dim;
    n1 = Dim.n1; m1 = Dim.m1;
    n2 = Dim.n2; m2 = Dim.m2;
    n3 = Dim.n3; m3 = Dim.m3;  
    
    phys = opt.phys;
    lattice_vec = phys.lattice.lattice_vector;
    [ Bodycenter_coef, Standard_coef,Edgecenter_coef_x,Edgecenter_coef_y, Edgecenter_coef_z,Facecenter_coef_x,Facecenter_coef_y,Facecenter_coef_z] = Grid_Generate_Mult_Point( Dim );
    
    switch phys.geometry.type
        case 'Isofunction'
            Standard_inner_idx  = Material_Locate_Isofunction( Standard_coef  , lattice_vec, phys);
            % filename = sprintf('./radius_data/radius_rate_square=%.3f.mat', radius_rate);
            % save(filename, 'Standard_inner_idx');
            Bodycenter_inner_idx = Material_Locate_Isofunction( Bodycenter_coef, lattice_vec, phys);
            Edgecenter_x_inner_index = Material_Locate_Isofunction( Edgecenter_coef_x, lattice_vec, phys);
            Edgecenter_y_inner_index = Material_Locate_Isofunction( Edgecenter_coef_y, lattice_vec, phys);
            Edgecenter_z_inner_index = Material_Locate_Isofunction( Edgecenter_coef_z, lattice_vec, phys);
            Facecenter_x_inner_index = Material_Locate_Isofunction( Facecenter_coef_x, lattice_vec, phys);
            Facecenter_y_inner_index = Material_Locate_Isofunction( Facecenter_coef_y, lattice_vec, phys);
            Facecenter_z_inner_index = Material_Locate_Isofunction( Facecenter_coef_z, lattice_vec, phys);
    end
   
    if strcmp(type, 'bianisotripic') % 都在顶点插值
        mtx_B.N = DistMaterialMatrix(m1, m2, m3, phys.permittivity_cov, Standard_inner_idx);
        mtx_B.M = DistMaterialMatrix(n1, n2, n3, phys.permeability_cov, Standard_inner_idx);
        %mtx_B.L = DistMaterialMatrix(n1, n2, n3, phys.magnetoelectric_zeta_cov, Standard_inner_idx);
        mtx_B.J = DistMaterialMatrix(n1, n2, n3, phys.magnetoelectric_Xi_cov, Standard_inner_idx);
    end



end

function A = DistMaterialMatrix(m1, m2, m3, par,idx)
% I_tilde = zeros(m1, m2, m3); I_tilde(idx) = 1;
n = m1*m2*m3;
I_tilde = sparse(idx, 1, 1, n, 1); 
for i = 1:3
    for j = 1:3       
    A{i,j} = spdiags( (par{2}(i,j) - par{1}(i,j))*I_tilde + par{1}(i,j), ...
                              0, n, n );
    end
end
end

