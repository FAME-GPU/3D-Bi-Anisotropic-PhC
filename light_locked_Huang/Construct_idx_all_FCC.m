function [edge_idx_in, edge_idx_out, Cubic_idx_in, Cubic_idx_out, face_idx_in, face_idx_out, vertex_idx_in, vertex_idx_out] = ...
    Construct_idx_all_FCC(n_1,n_2,n_3,edge_len,rd_sphere,rd_cylinder, delta,shape)
 
                
n        = n_1 * n_2 * n_3; 
xyz      = zeros(n,3);

xx       = (0.5:n_1-0.5)' * delta(1);
xyz(:,1) = kron(ones(n_2*n_3,1),xx);

yy       = (0:n_2-1)' * delta(2);
yy1      = kron(yy,ones(n_1,1));
xyz(:,2) = kron(ones(n_3,1),yy1);

zz       = (0:n_3-1)' * delta(3);
xyz(:,3) = kron(zz,ones(n_1*n_2,1));

switch shape
    case 'cylinder'
        InOut = locate_cyld(xyz, rd_cylinder, rd_sphere, edge_len);
    case 'ellipse'
        InOut = locate_ellipse(xyz, rd_sphere, rd_cylinder, edge_len);
end

edge_idx_in.X  = find(InOut == 1);
edge_idx_out.X = find(InOut == 0);

xx       = (0:n_1-1)' * delta(1);
xyz(:,1) = kron(ones(n_2*n_3,1),xx);

yy       = (0.5:n_2-0.5)' * delta(2);
yy1      = kron(yy,ones(n_1,1));
xyz(:,2) = kron(ones(n_3,1),yy1);

zz       = (0:n_3-1)' * delta(3);
xyz(:,3) = kron(zz,ones(n_1*n_2,1));

switch shape
    case 'cylinder'
        InOut = locate_cyld(xyz, rd_cylinder, rd_sphere, edge_len);
    case 'ellipse'
        InOut = locate_ellipse(xyz, rd_sphere, rd_cylinder, edge_len);
end

edge_idx_in.Y    = find(InOut == 1);
edge_idx_out.Y   = find(InOut == 0); 

xx       = (0:n_1-1)' * delta(1);
xyz(:,1) = kron(ones(n_2*n_3,1),xx);

yy       = (0:n_2-1)' * delta(2);
yy1      = kron(yy,ones(n_1,1));
xyz(:,2) = kron(ones(n_3,1),yy1);

zz       = (0.5:n_3-0.5)' * delta(3);
xyz(:,3) = kron(zz,ones(n_1*n_2,1));

switch shape
    case 'cylinder'
        InOut = locate_cyld(xyz, rd_cylinder, rd_sphere, edge_len);
    case 'ellipse'
        InOut = locate_ellipse(xyz, rd_sphere, rd_cylinder, edge_len);
end

edge_idx_in.Z    = find(InOut == 1);
edge_idx_out.Z   = find(InOut == 0);

%
xx       = (0.5:n_1-0.5)' * delta(1);
xyz(:,1) = kron(ones(n_2*n_3,1),xx);

yy       = (0.5:n_2-0.5)' * delta(2);
yy1      = kron(yy,ones(n_1,1));
xyz(:,2) = kron(ones(n_3,1),yy1);

zz       = (0.5:n_3-0.5)' * delta(3);
xyz(:,3) = kron(zz,ones(n_1*n_2,1));

switch shape
    case 'cylinder'
        InOut = locate_cyld(xyz, rd_cylinder, rd_sphere, edge_len);
    case 'ellipse'
        InOut = locate_ellipse(xyz, rd_sphere, rd_cylinder, edge_len);
end

Cubic_idx_in       = find(InOut == 1);
Cubic_idx_out      = find(InOut == 0);

%
% Face
%
%
% (I)
%
xx       = (0:n_1-1)' * delta(1);
xyz(:,1) = kron(ones(n_2*n_3,1),xx);

yy       = (0.5:n_2-0.5)' * delta(2);
yy1      = kron(yy,ones(n_1,1));
xyz(:,2) = kron(ones(n_3,1),yy1);

zz       = (0.5:n_3-0.5)' * delta(3);
xyz(:,3) = kron(zz,ones(n_1*n_2,1));

switch shape
    case 'cylinder'
        InOut = locate_cyld(xyz, rd_cylinder, rd_sphere, edge_len);
    case 'ellipse'
        InOut = locate_ellipse(xyz, rd_sphere, rd_cylinder, edge_len);
end

face_idx_in.X  = find(InOut == 1);
face_idx_out.X = find(InOut == 0); 

%
% (II)
%
xx       = (0.5:n_1-0.5)' * delta(1);
xyz(:,1) = kron(ones(n_2*n_3,1),xx);

yy       = (0:n_2-1)' * delta(2);
yy1      = kron(yy,ones(n_1,1));
xyz(:,2) = kron(ones(n_3,1),yy1);

zz       = (0.5:n_3-0.5)' * delta(3);
xyz(:,3) = kron(zz,ones(n_1*n_2,1));

switch shape
    case 'cylinder'
        InOut = locate_cyld(xyz, rd_cylinder, rd_sphere, edge_len);
    case 'ellipse'
        InOut = locate_ellipse(xyz, rd_sphere, rd_cylinder, edge_len);
end

face_idx_in.Y    = find(InOut == 1);
face_idx_out.Y   = find(InOut == 0); 

%
% (III)
%
xx       = (0.5:n_1-0.5)' * delta(1);
xyz(:,1) = kron(ones(n_2*n_3,1),xx);

yy       = (0.5:n_2-0.5)' * delta(2);
yy1      = kron(yy,ones(n_1,1));
xyz(:,2) = kron(ones(n_3,1),yy1);

zz       = (0:n_3-1)' * delta(3);
xyz(:,3) = kron(zz,ones(n_1*n_2,1));

switch shape
    case 'cylinder'
        InOut = locate_cyld(xyz, rd_cylinder, rd_sphere, edge_len);
    case 'ellipse'
        InOut = locate_ellipse(xyz, rd_sphere, rd_cylinder, edge_len);
end

face_idx_in.Z    = find(InOut == 1);
face_idx_out.Z   = find(InOut == 0);

%
xx       = (0:n_1-1)' * delta(1);
xyz(:,1) = kron(ones(n_2*n_3,1),xx);

yy       = (0:n_2-1)' * delta(2);
yy1      = kron(yy,ones(n_1,1));
xyz(:,2) = kron(ones(n_3,1),yy1);

zz       = (0:n_3-1)' * delta(3);
xyz(:,3) = kron(zz,ones(n_1*n_2,1));

switch shape
    case 'cylinder'
        InOut = locate_cyld(xyz, rd_cylinder, rd_sphere, edge_len);
    case 'ellipse'
        InOut = locate_ellipse(xyz, rd_sphere, rd_cylinder, edge_len);
end

vertex_idx_in    = find(InOut == 1);
vertex_idx_out   = find(InOut == 0);

end

% =========================================================================
function InOut = locate_ellipse(xyz, rd_sphere, rd_ellipse, edge_len)
        
    shift   = [ 0.0, 0.0, 0.0 ];
    
    S1 = edge_len * [0 ,0 ,0] + shift;                      % (0, 0, 0)
    S2 = edge_len * [1 ,0 ,0] + shift;                      % a_1
    S3 = edge_len * [1/2 ,sqrt(3)/2 ,0        ] + shift;    % a_2
    S4 = edge_len * [1/2 ,sqrt(3)/6 ,sqrt(6)/3] + shift;    % a_3
    S5 = edge_len * [1/2 ,sqrt(3)/6 ,sqrt(6)/12]+ shift;    % ¤¤¶¡¨ºÁû
    
    n          = size(xyz,1);
    InOut      = zeros(n,1);
    
    xyz_S5     = xyz-kron(S5,ones(n,1));
    nrm_xyz_S5 = sum(xyz_S5.*xyz_S5,2);
    nrm_xyz_S5 = nrm_xyz_S5.^(0.5);    

%     if     (norm(point-S1)+norm(point-S5)) <= 2*sqrt(rd_ellipse^2+(norm(S1-S5)/2)^2)
%         in_out = 'I';
    xyz_S5     = xyz-kron(S1,ones(n,1));
    nrm_xyz_S1 = sum(xyz_S5.*xyz_S5,2);
    nrm_xyz_S1 = nrm_xyz_S1.^(0.5)+nrm_xyz_S5;
    idx        = find( nrm_xyz_S1 <= 2*sqrt(rd_ellipse^2+(norm(S1-S5)/2)^2) );
    InOut(idx) = ones(length(idx),1);
        
%     elseif (norm(point-S2)+norm(point-S5)) <= 2*sqrt(rd_ellipse^2+(norm(S2-S5)/2)^2)
%         in_out = 'I';
    xyz_S5     = xyz-kron(S2,ones(n,1));
    nrm_xyz_S1 = sum(xyz_S5.*xyz_S5,2);
    nrm_xyz_S1 = nrm_xyz_S1.^(0.5)+nrm_xyz_S5;
    idx        = find( nrm_xyz_S1 <= 2*sqrt(rd_ellipse^2+(norm(S2-S5)/2)^2) );
    InOut(idx) = ones(length(idx),1);
    
%     elseif (norm(point-S3)+norm(point-S5)) <= 2*sqrt(rd_ellipse^2+(norm(S3-S5)/2)^2)
%         in_out = 'I';
    xyz_S5     = xyz-kron(S3,ones(n,1));
    nrm_xyz_S1 = sum(xyz_S5.*xyz_S5,2);
    nrm_xyz_S1 = nrm_xyz_S1.^(0.5)+nrm_xyz_S5;
    idx        = find( nrm_xyz_S1 <= 2*sqrt(rd_ellipse^2+(norm(S3-S5)/2)^2) );
    InOut(idx) = ones(length(idx),1);
    
%     elseif (norm(point-S4)+norm(point-S5)) <= 2*sqrt(rd_ellipse^2+(norm(S4-S5)/2)^2)
%         in_out = 'I';
    xyz_S5     = xyz-kron(S4,ones(n,1));
    nrm_xyz_S1 = sum(xyz_S5.*xyz_S5,2);
    nrm_xyz_S1 = nrm_xyz_S1.^(0.5)+nrm_xyz_S5;
    idx        = find( nrm_xyz_S1 <= 2*sqrt(rd_ellipse^2+(norm(S4-S5)/2)^2) );
    InOut(idx) = ones(length(idx),1);
    
%     end

%     if in_out ~= 'I' && norm(point-S1) <= rd_sphere   % S1
%         in_out = 'I';
%
xyz_new    = xyz - kron(S1,ones(n,1));
D1         = sum(xyz_new.*xyz_new,2);
idx        = find( D1 <= rd_sphere^2 );
InOut(idx) = ones(length(idx),1);

%     elseif norm(point-S2) <= rd_sphere                % S2
%         in_out = 'I';
%
xyz_new    = xyz - kron(S2,ones(n,1));
D1         = sum(xyz_new.*xyz_new,2);
idx        = find( D1 <= rd_sphere^2 );
InOut(idx) = ones(length(idx),1);

%     elseif norm(point-S3) <= rd_sphere                % S3
%         in_out = 'I';
%
xyz_new    = xyz - kron(S3,ones(n,1));
D1         = sum(xyz_new.*xyz_new,2);
idx        = find( D1 <= rd_sphere^2 );
InOut(idx) = ones(length(idx),1);

%     elseif norm(point-S4) <= rd_sphere                % S4
%         in_out = 'I';
%
xyz_new    = xyz - kron(S4,ones(n,1));
D1         = sum(xyz_new.*xyz_new,2);
idx        = find( D1 <= rd_sphere^2 );
InOut(idx) = ones(length(idx),1);

%     elseif norm(point-S5) <= rd_sphere                % S5
%         in_out = 'I';
%
xyz_new    = xyz - kron(S5,ones(n,1));
D1         = sum(xyz_new.*xyz_new,2);
idx        = find( D1 <= rd_sphere^2 );
InOut(idx) = ones(length(idx),1);

end

% =========================================================================
function InOut = locate_cyld(xyz, rd_cylinder, rd_sphere, edge_len)

shift     = [ 0.0, 0.0, 0.0 ]; %-edge_len * [1/2 ,sqrt(3)/6 ,sqrt(6)/12];  %[ 0.0, 0.0, 0.0 ];

S1        = edge_len * [0 ,0 ,0] + shift;                      % (0, 0, 0)
S2        = edge_len * [1 ,0 ,0] + shift;                      % a_1
S3        = edge_len * [1/2 ,sqrt(3)/2 ,0        ] + shift;    % a_2
S4        = edge_len * [1/2 ,sqrt(3)/6 ,sqrt(6)/3] + shift;    % a_3
S5        = edge_len * [1/2 ,sqrt(3)/6 ,sqrt(6)/12]+ shift;    % ????????

distS15   = norm(S1-S5);
distS25   = norm(S2-S5);
distS35   = norm(S3-S5);
distS45   = norm(S4-S5);

tmpS15    = (S1-S5)/distS15;
tmpS25    = (S2-S5)/distS25;
tmpS35    = (S3-S5)/distS35;
tmpS45    = (S4-S5)/distS45;

uuT15     = tmpS15.' * tmpS15;
uuT25     = tmpS25.' * tmpS25;
uuT35     = tmpS35.' * tmpS35;
uuT45     = tmpS45.' * tmpS45;

% projToS15 = (point-S5)*uuT15.' + S5;
% projToS25 = (point-S5)*uuT25.' + S5;
% projToS35 = (point-S5)*uuT35.' + S5;
% projToS45 = (point-S5)*uuT45.' + S5;
 
n          = size(xyz,1);
xyz_new    = xyz-kron(S5,ones(n,1));
InOut      = zeros(n,1);

%     if     norm(point-projToS15) <= rd_cylinder && norm(projToS15-S1) <= distS15 && norm(projToS15-S5) <= distS15
%         in_out = 'I';
% 
projToS15  = xyz_new*uuT15.' + kron(S5,ones(n,1));
tmp1       = xyz - projToS15;
D1         = sum(tmp1.*tmp1,2);
tmp1       = projToS15 - kron(S1,ones(n,1));
D2         = sum(tmp1.*tmp1,2);
tmp1       = projToS15 - kron(S5,ones(n,1));
D3         = sum(tmp1.*tmp1,2);
idx        = find( D1 <= rd_cylinder^2 & D2 <= distS15^2 & D3 <= distS15^2 );
InOut(idx) = ones(length(idx),1);

%     elseif norm(point-projToS25) <= rd_cylinder && norm(projToS25-S2) <= distS25 && norm(projToS25-S5) <= distS25
%         in_out = 'I';
% 
projToS15  = xyz_new*uuT25.' + kron(S5,ones(n,1));
tmp1       = xyz - projToS15;
D1         = sum(tmp1.*tmp1,2);
tmp1       = projToS15 - kron(S2,ones(n,1));
D2         = sum(tmp1.*tmp1,2);
tmp1       = projToS15 - kron(S5,ones(n,1));
D3         = sum(tmp1.*tmp1,2);
idx        = find( D1 <= rd_cylinder^2 & D2 <= distS25^2 & D3 <= distS25^2 );
InOut(idx) = ones(length(idx),1);

%     elseif norm(point-projToS35) <= rd_cylinder && norm(projToS35-S3) <= distS35 && norm(projToS35-S5) <= distS35
%         in_out = 'I';
%         
projToS15  = xyz_new*uuT35.' + kron(S5,ones(n,1));
tmp1       = xyz - projToS15;
D1         = sum(tmp1.*tmp1,2);
tmp1       = projToS15 - kron(S3,ones(n,1));
D2         = sum(tmp1.*tmp1,2);
tmp1       = projToS15 - kron(S5,ones(n,1));
D3         = sum(tmp1.*tmp1,2);
idx        = find( D1 <= rd_cylinder^2 & D2 <= distS35^2 & D3 <= distS35^2 );
InOut(idx) = ones(length(idx),1);

%     elseif norm(point-projToS45) <= rd_cylinder && norm(projToS45-S4) <= distS45 && norm(projToS45-S5) <= distS45
%         in_out = 'I';  
%         
projToS15  = xyz_new*uuT45.' + kron(S5,ones(n,1));
tmp1       = xyz - projToS15;
D1         = sum(tmp1.*tmp1,2);
tmp1       = projToS15 - kron(S4,ones(n,1));
D2         = sum(tmp1.*tmp1,2);
tmp1       = projToS15 - kron(S5,ones(n,1));
D3         = sum(tmp1.*tmp1,2);
idx        = find( D1 <= rd_cylinder^2 & D2 <= distS45^2 & D3 <= distS45^2 );
InOut(idx) = ones(length(idx),1);
% 
%     end
 
%     if in_out ~= 'I' && norm(point-S1) <= rd_sphere   % S1
%         in_out = 'I';
%
xyz_new    = xyz - kron(S1,ones(n,1));
D1         = sum(xyz_new.*xyz_new,2);
idx        = find( D1 <= rd_sphere^2 );
InOut(idx) = ones(length(idx),1);

%     elseif norm(point-S2) <= rd_sphere                % S2
%         in_out = 'I';
%
xyz_new    = xyz - kron(S2,ones(n,1));
D1         = sum(xyz_new.*xyz_new,2);
idx        = find( D1 <= rd_sphere^2 );
InOut(idx) = ones(length(idx),1);

%     elseif norm(point-S3) <= rd_sphere                % S3
%         in_out = 'I';
%
xyz_new    = xyz - kron(S3,ones(n,1));
D1         = sum(xyz_new.*xyz_new,2);
idx        = find( D1 <= rd_sphere^2 );
InOut(idx) = ones(length(idx),1);

%     elseif norm(point-S4) <= rd_sphere                % S4
%         in_out = 'I';
%
xyz_new    = xyz - kron(S4,ones(n,1));
D1         = sum(xyz_new.*xyz_new,2);
idx        = find( D1 <= rd_sphere^2 );
InOut(idx) = ones(length(idx),1);

%     elseif norm(point-S5) <= rd_sphere                % S5
%         in_out = 'I';
%
xyz_new    = xyz - kron(S5,ones(n,1));
D1         = sum(xyz_new.*xyz_new,2);
idx        = find( D1 <= rd_sphere^2 );
InOut(idx) = ones(length(idx),1);

end
