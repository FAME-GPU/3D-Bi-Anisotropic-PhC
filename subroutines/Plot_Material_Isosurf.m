function Plot_Material_Isosurf(opt)
% plot settings
figure(11); cla
hold on
color_map{1} = [1,0,0];
color_map{2} = [0,0,1];
geometry = opt.phys.geometry;

Dim = opt.comp.Dim;
n1 = Dim.n1;  n2 = Dim.n2;  n3 = Dim.n3;

phys = opt.phys;
lattice_vec = phys.lattice.lattice_vector;

[ Bodycenter_coef, ~ ] = Grid_Generate( Dim );

supercell_lattice_vec = lattice_vec * diag(geometry.supercell_num);
Point_set = Bodycenter_coef*supercell_lattice_vec';
unitcell_coef = Point_set/lattice_vec';

X = permute(reshape(Point_set(:,1),n1,n2,n3),[2,1,3]);
Y = permute(reshape(Point_set(:,2),n1,n2,n3),[2,1,3]);
Z = permute(reshape(Point_set(:,3),n1,n2,n3),[2,1,3]);

for i = 1:length(geometry.isofun)
    f = geometry.isofun{i};
    v = geometry.isoval(i);
    V{i} = f(Point_set,unitcell_coef);
    if isfield(geometry,'mask')
        mask = geometry.mask{i};
        mask_idx = find(mask(Point_set,unitcell_coef) == 1);
        V{i}(mask_idx) = 0;
    end
    V{i} = permute(reshape(V{i},n1,n2,n3),[2,1,3]);
    p_surf{i} = patch( isosurface(X,Y,Z,V{i},v,'r') );
    p_caps{i} = patch( isocaps(X,Y,Z,V{i},v) );
    p_surf{i}.FaceColor = color_map{i};
    p_caps{i}.FaceColor = 'blue';
    p_surf{i}.EdgeColor = 'none';
    p_caps{i}.EdgeColor = 'none';
    % p{i}.FaceAlpha = 0.5;
end
camlight right
rotate3d on
axis equal
axis tight
axis off
view(90,0)

%% Plot edge of each unit cell
n1 = geometry.supercell_num(1)+1;
n2 = geometry.supercell_num(2)+1;
n3 = geometry.supercell_num(3)+1;
x = (0:n1-1)';
y = (0:n2-1)';% y = (0:n2-1)';
z = (0:n3-1)';       
X = kron( ones(n3,1), kron(ones(n2,1),         x) );
Y = kron( ones(n3,1), kron(         y,ones(n1,1)) );
Z = kron(          z, kron(ones(n2,1),ones(n1,1)) );
Ax = diag(ones(n1-1,1),1) + diag(ones(n1,1),0);
Ay = diag(ones(n2-1,1),1) + diag(ones(n2,1),0);
Az = diag(ones(n3-1,1),1) + diag(ones(n3,1),0);
A  = kron(eye(n3),kron(eye(n2),Ax'*Ax)) + kron(eye(n3),kron(Ay'*Ay,eye(n1))) + kron(Az'*Az,kron(eye(n2),eye(n1))) ;
A  = triu(A);
[pt1,pt2] = ind2sub([n1*n2*n3,n1*n2*n3],find(A==1));
edge = [pt1,pt2];

Node_set = [X,Y,Z]*lattice_vec';

for i = 1:size(edge,1)
    plot3(Node_set(edge(i,:),1),Node_set(edge(i,:),2),Node_set(edge(i,:),3),'k-');
end


% %% Compute Rotation matrix s.t a2 -> y-axis
% r = vrrotvec(lattice_vec(:,2),norm(lattice_vec(:,2))*[0,1,0]');
% R = vrrotvec2mat(r);
% RotPoint_set = (R*Point_set')';
% %% Generate standard 3D grid
% Pad_num = 3;
% x_lim = [ min(RotPoint_set(:,1)), max(RotPoint_set(:,1)) ] + [-Pad_num,Pad_num]/n;
% y_lim = [ min(RotPoint_set(:,2)), max(RotPoint_set(:,2)) ] + [-Pad_num,Pad_num]/n;
% z_lim = [ min(RotPoint_set(:,3)), max(RotPoint_set(:,3)) ] + [-Pad_num,Pad_num]/n;
% node_num = [n,n,n].*geometry.supercell_num;
% x = linspace(x_lim(1),x_lim(2),node_num(1));
% y = linspace(y_lim(1),y_lim(2),node_num(2));
% z = linspace(z_lim(1),z_lim(2),node_num(3));
% [X,Y,Z] = meshgrid(x,y,z);
% Node_set = [X(:),Y(:),Z(:)];
% Pad_idx = find( X(:)<min(RotPoint_set(:,1)) | X(:)>max(RotPoint_set(:,1)) |...
%                 Y(:)<min(RotPoint_set(:,2)) | Y(:)>max(RotPoint_set(:,2)) |...
%                 Z(:)<min(RotPoint_set(:,3)) | Z(:)>max(RotPoint_set(:,3)) );
% %% Cancel points below span(a1,a3) and above span(a1+k*a2,a3+k*ax) (k2:=geometry.supercell_num(2))
% a1 = R*lattice_vec(:,1);
% a2 = R*lattice_vec(:,2);
% a3 = R*lattice_vec(:,3);
% % Below span(a1,a3)
% normvec13 = cross(a1,a3);
% Pad_idx = union( Pad_idx, find( Node_set*normvec13 >0 ) );
% % Above span(a1+k*a2,a3+k*ax)
% tmp_o = geometry.supercell_num(2)*a2;
% Pad_idx = union( Pad_idx, find( (Node_set-tmp_o')*normvec13 <0 ) );
% % Below span(a1,a2)
% normvec12 = cross(a1,a2);
% Pad_idx = union( Pad_idx, find( Node_set*normvec12 <0 ) );
% % Above span(a1,a2)
% tmp_o = geometry.supercell_num(1)*a1 + geometry.supercell_num(3)*a3;
% Pad_idx = union( Pad_idx, find(( Node_set-tmp_o')*normvec12 > 0 ) );
% % Below span(a2,a3)
% normvec23 = cross(a2,a3);
% Pad_idx = union( Pad_idx, find( Node_set*normvec23 <0 ) );
% % Above span(a2,a3)
% tmp_o = geometry.supercell_num(1)*a1 + geometry.supercell_num(3)*a3;
% Pad_idx = union( Pad_idx, find(( Node_set-tmp_o')*normvec23 > 0 ) );
% %% Construct isosurface
% for i = 1:length(geometry.isofun)
%     f = geometry.isofun{i};
%     v = geometry.isoval(i);
%     RotNode_set = (R'*Node_set')';
%     Node_coef = RotNode_set/lattice_vec';
%     V{i} = f(RotNode_set,Node_coef);
%     if isfield(geometry,'mask')
%         mask = geometry.mask{i};
%         mask_idx = find(mask(RotNode_set, Node_coef) == 1);
%         V{i}(mask_idx) = nan;
%     end
% 
%     V{i}(Pad_idx) = nan;
%     V{i} = reshape(V{i},node_num([2,1,3]));
% 
% 
%     p{i} = patch( isosurface(X,Y,Z,V{i},v,'r') );
%     % p{i} = patch( isocaps(X,Y,Z,V{i},v) );
%     p{i}.FaceColor = color_map{i};
%     % p{i}.FaceColor = 'interp';
%     p{i}.EdgeColor = 'none';
%     % p{i}.FaceAlpha = 0.5;
% end
% camlight right
% rotate3d on
% axis equal
% axis tight
% axis off
% view(90,0)
% %% Plot lattice vectors
% hax = gca;
% O = [0;0;0] + tau*a2;
% stemWidth = 0.01*norm(lattice_vec(:,1));
% tipWidth  = 0.025*norm(lattice_vec(:,1));
% mArrow3(hax,O',(R*lattice_vec(:,1))'+tau*a2','color',  'red','stemWidth',stemWidth,'tipWidth',tipWidth,'facealpha',0.9);
% mArrow3(hax,O',(R*lattice_vec(:,2))'+tau*a2','color','green','stemWidth',stemWidth,'tipWidth',tipWidth,'facealpha',0.9);
% mArrow3(hax,O',(R*lattice_vec(:,3))'+tau*a2','color', 'blue','stemWidth',stemWidth,'tipWidth',tipWidth,'facealpha',0.9);
% % plot3(RotPoint_set(Point_idx,1),RotPoint_set(Point_idx,2),RotPoint_set(Point_idx,3),'b.');
% %% Plot edge of each unit cell
% n1 = geometry.supercell_num(1)+1;
% n2 = geometry.supercell_num(2)+2;
% n3 = geometry.supercell_num(3)+1;
% x = (0:n1-1)';
% y = [0,(0:n2-3)+tau, n2-2]';% y = (0:n2-1)';
% z = (0:n3-1)';       
% X = kron( ones(n3,1), kron(ones(n2,1),         x) );
% Y = kron( ones(n3,1), kron(         y,ones(n1,1)) );
% Z = kron(          z, kron(ones(n2,1),ones(n1,1)) );
% Ax = diag(ones(n1-1,1),1) + diag(ones(n1,1),0);
% Ay = diag(ones(n2-1,1),1) + diag(ones(n2,1),0);
% Az = diag(ones(n3-1,1),1) + diag(ones(n3,1),0);
% A  = kron(eye(n3),kron(eye(n2),Ax'*Ax)) + kron(eye(n3),kron(Ay'*Ay,eye(n1))) + kron(Az'*Az,kron(eye(n2),eye(n1))) ;
% A  = triu(A);
% [pt1,pt2] = ind2sub([n1*n2*n3,n1*n2*n3],find(A==1));
% edge = [pt1,pt2];
% 
% Node_set = [X,Y,Z]*[a1,a2,a3]';
% 
% for i = 1:size(edge,1)
%     plot3(Node_set(edge(i,:),1),Node_set(edge(i,:),2),Node_set(edge(i,:),3),'k-');
% end
% axis equal
% %% Plot termination plane
% O = ceil(geometry.supercell_num(2)/2)*a2 + [0;0;0];
% PT = [O,O+a1,O+a1+a3,O+a3,O]';
% fill3(PT(:,1),PT(:,2),PT(:,3),'r','FaceAlpha',0.5);
% 
% %% Plot nodes in material
% plot3(RotPoint_set(Point_idx,1),RotPoint_set(Point_idx,2),RotPoint_set(Point_idx,3),'k.');
