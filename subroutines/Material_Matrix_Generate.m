function [ mtx_B ] = Material_Matrix_Generate( opt )
    Dim = opt.comp.Dim;
    n1 = Dim.n1; m1 = Dim.m1;
    n2 = Dim.n2; m2 = Dim.m2;
    n3 = Dim.n3; m3 = Dim.m3;  
    
    phys = opt.phys;
    lattice_vec = phys.lattice.lattice_vector;

    [ Bodycenter_coef, Standard_coef ] = Grid_Generate( Dim );
    
    switch phys.geometry.type
        case 'Isofunction'
            Standard_inner_idx   = Material_Locate_Isofunction( Standard_coef  , lattice_vec, phys.geometry);
            Bodycenter_inner_idx = Material_Locate_Isofunction( Bodycenter_coef, lattice_vec, phys.geometry);
    end
    
    % mtx_B.N = DistMaterialMatrix(m1, m2, m3, phys.inv_permittivity_cov,   Standard_inner_idx);
    % mtx_B.M = DistMaterialMatrix(n1, n2, n3, phys.inv_permeability_cov, Bodycenter_inner_idx); 
    mtx_B.Standard_inner_idx   = Standard_inner_idx;
    mtx_B.Bodycenter_inner_idx = Bodycenter_inner_idx;
    
    %% 手性
    mtx_B.N = DistMaterialMatrix(m1, m2, m3, phys.permittivity_cov,   Standard_inner_idx);  % 
    mtx_B.M = DistMaterialMatrix(n1, n2, n3, phys.permeability_cov, Bodycenter_inner_idx); 
    mtx_B.L = DistMaterialMatrix(m1, m2, m3, phys.magnetoelectric_zeta_cov,   Standard_inner_idx);
    mtx_B.J = DistMaterialMatrix(m1, m2, m3, phys.magnetoelectric_Xi_cov,   Standard_inner_idx);
    mtx_B.M_inv = DistMaterialMatrix(m1, m2, m3, phys.permeability_cov_inv,   Standard_inner_idx);
    mtx_B.Phi_inv = DistMaterialMatrix(m1, m2, m3, phys.Phi_inv,   Standard_inner_idx);
end

function A = DistMaterialMatrix(n1, n2, n3, par, inner_idx)
    I_tilde = zeros(n1,n2,n3); I_tilde(inner_idx) = 1;
    n = n1*n2*n3;
    for i = 1:3
        for j = 1:3
            A{i,j} = spdiags( (par{2}(i,j) - par{1}(i,j))*I_tilde(:) + par{1}(i,j), ...
                              0, n, n );
        end
    end
end


function plotMaterial(Dim,edge_len,lattice_vec,Bodycenter_point_set,Bodycenter_inner_idx)
    O = -0.5*edge_len';

    a1 = O + lattice_vec(:,1);
    a2 = O + lattice_vec(:,2);
    a3 = O + lattice_vec(:,3);

    hax = gca;
    mArrow3(hax,O',a1','color','red','stemWidth',0.15,'facealpha',0.5); hold on
    mArrow3(hax,O',a2','color','green','stemWidth',0.15,'facealpha',0.5);
    mArrow3(hax,O',a3','color','blue','stemWidth',0.15,'facealpha',0.5);

    plot3(Bodycenter_point_set(:,1),...
          Bodycenter_point_set(:,2),...
          Bodycenter_point_set(:,3), 'b.');
    axis tight
    axis equal
%     axis([0,edge_len(1),0,edge_len(2),0,edge_len(3)])

    plot3(Bodycenter_point_set(Bodycenter_inner_idx,1),...
          Bodycenter_point_set(Bodycenter_inner_idx,2),...
          Bodycenter_point_set(Bodycenter_inner_idx,3), 'ro');
%     axis tight
%     axis equal
%     axis([0,edge_len(1),0,edge_len(2),0,edge_len(3)])
%     subplot(1,3,3)
% 
%     figure(10);
%     hax = gca;
%     mArrow3(hax,O',a1','color','red','stemWidth',0.15,'facealpha',0.5); hold on
%     mArrow3(hax,O',a2','color','green','stemWidth',0.15,'facealpha',0.5);
%     mArrow3(hax,O',a3','color','blue','stemWidth',0.15,'facealpha',0.5);
% 
%     x = reshape(Bodycenter_point_set(:,1),Dim.row_1D.x,Dim.row_1D.y,Dim.row_1D.z);
%     y = reshape(Bodycenter_point_set(:,2),Dim.row_1D.x,Dim.row_1D.y,Dim.row_1D.z);
%     z = reshape(Bodycenter_point_set(:,3),Dim.row_1D.x,Dim.row_1D.y,Dim.row_1D.z);
%     v = zeros(size(x));
%     v(Standard_inner_idx) = 1;
% 
%     X = zeros(size(x)+2); X(2:end-1,2:end-1,2:end-1) = x; X(1,:,:) = min(x(:)) - mesh_len(1); X(end,:,:) = max(x(:)) + mesh_len(1);
%     Y = zeros(size(y)+2); Y(2:end-1,2:end-1,2:end-1) = y; Y(1,:,:) = min(y(:)) - mesh_len(2); Y(end,:,:) = max(y(:)) + mesh_len(2);
%     Z = zeros(size(z)+2); Z(2:end-1,2:end-1,2:end-1) = z; Z(1,:,:) = min(z(:)) - mesh_len(3); Z(end,:,:) = max(z(:)) + mesh_len(3);
%     V = zeros(size(v)+2); V(2:end-1,2:end-1,2:end-1) = v;
% 
%     p=patch(isosurface(X,Y,Z,V,0.99));
% %         p=patch(isosurface(x,y,z,v,0.99));
% %         isonormals(x,y,z,v,p)
%     p.FaceColor = [0, 229, 255]/255;
%     p.EdgeColor = 'none';
%     daspect([1 1 1])
%     view(3); 
%     axis([0,edge_len(1),0,edge_len(2),0,edge_len(3)]);
% %         axis tight
%     camlight 
%     lighting gouraud 

end

function plotMesh(edge_len,grid_num)
    fontSize=1;
    faceAlpha1=0;

%     Nx = 4; Ny = 4; Nz = 4;

    StdboxSize    = edge_len;
    StdboxCellNum = grid_num;

    [meshStruct]=hexMeshBox(StdboxSize,StdboxCellNum);

    E  = meshStruct.E;
    V  = meshStruct.V;
    F  = meshStruct.F;
    Fb = meshStruct.Fb;
    faceBoundaryMarker = meshStruct.faceBoundaryMarker;

    cFigure;
    hold on;

    gpatch(Fb,V,faceBoundaryMarker,'k',faceAlpha1);
    % patchNormPlot(Fb,V);

    axisGeom(gca,fontSize);
    colormap(gjet(6)); icolorbar;
    drawnow;
end