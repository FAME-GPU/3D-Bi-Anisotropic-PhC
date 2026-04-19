function [mtx_parameter] = Construct_parameter_mtx(n, face_edge_idx_in, vertex_Cubic_idx_in, anisotropic_mtx, value_air)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
i_idx = zeros(36*n,1);
j_idx = zeros(36*n,1);
val   = zeros(36*n,1);
%
% block (1,1)
%

tmp_x                       = value_air(1,1) * ones(n,1);
tmp_x(face_edge_idx_in.X,1) = anisotropic_mtx(1,1);
i_idx(1:n,1)                = (1:n)';
j_idx(1:n,1)                = (1:n)';
val(1:n,1)                  = tmp_x;
leng                        = n;

tmp_x                       = value_air(2,2) * ones(n,1);
tmp_x(face_edge_idx_in.Y,1) = anisotropic_mtx(2,2);
i_idx(leng+1:leng+n,1)      = (n+1:2*n)';
j_idx(leng+1:leng+n,1)      = (n+1:2*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(3,3) * ones(n,1);
tmp_x(face_edge_idx_in.Z,1) = anisotropic_mtx(3,3);
i_idx(leng+1:leng+n,1)      = (2*n+1:3*n)';
j_idx(leng+1:leng+n,1)      = (2*n+1:3*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

%
% block (1,2)
%
tmp_x                       = value_air(1,2) * ones(n,1);
tmp_x(face_edge_idx_in.X,1) = anisotropic_mtx(1,2);
i_idx(leng+1:leng+n,1)      = (1:n)';
j_idx(leng+1:leng+n,1)      = (4*n+1:5*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(2,1) * ones(n,1);
tmp_x(face_edge_idx_in.Y,1) = anisotropic_mtx(2,1);
i_idx(leng+1:leng+n,1)      = (n+1:2*n)';
j_idx(leng+1:leng+n,1)      = (3*n+1:4*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

%
% block (1,3)
%
tmp_x                       = value_air(1,3) * ones(n,1);
tmp_x(face_edge_idx_in.X,1) = anisotropic_mtx(1,3);
i_idx(leng+1:leng+n,1)      = (1:n)';
j_idx(leng+1:leng+n,1)      = (8*n+1:9*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(3,1) * ones(n,1);
tmp_x(face_edge_idx_in.Z,1) = anisotropic_mtx(3,1);
i_idx(leng+1:leng+n,1)      = (2*n+1:3*n)';
j_idx(leng+1:leng+n,1)      = (6*n+1:7*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

%
% block (1,4)
%
tmp_x                       = value_air(2,3) * ones(n,1);
tmp_x(face_edge_idx_in.Y,1) = anisotropic_mtx(2,3);
i_idx(leng+1:leng+n,1)      = (n+1:2*n)';
j_idx(leng+1:leng+n,1)      = (11*n+1:12*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(3,2) * ones(n,1);
tmp_x(face_edge_idx_in.Z,1) = anisotropic_mtx(3,2);
i_idx(leng+1:leng+n,1)      = (2*n+1:3*n)';
j_idx(leng+1:leng+n,1)      = (10*n+1:11*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

% ======================================================
%
% block (2,1)
%
tmp_x                       = value_air(1,2) * ones(n,1);
tmp_x(face_edge_idx_in.Y,1) = anisotropic_mtx(1,2);
i_idx(leng+1:leng+n,1)      = (3*n+1:4*n)';
j_idx(leng+1:leng+n,1)      = (n+1:2*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(2,1) * ones(n,1);
tmp_x(face_edge_idx_in.X,1) = anisotropic_mtx(2,1);
i_idx(leng+1:leng+n,1)      = (4*n+1:5*n)';
j_idx(leng+1:leng+n,1)      = (1:n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

%
% block (2,2)
%
tmp_x                       = value_air(1,1) * ones(n,1);
tmp_x(face_edge_idx_in.Y,1) = anisotropic_mtx(1,1);
i_idx(leng+1:leng+n,1)      = (3*n+1:4*n)';
j_idx(leng+1:leng+n,1)      = (3*n+1:4*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(2,2) * ones(n,1);
tmp_x(face_edge_idx_in.X,1) = anisotropic_mtx(2,2);
i_idx(leng+1:leng+n,1)      = (4*n+1:5*n)';
j_idx(leng+1:leng+n,1)      = (4*n+1:5*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                        = value_air(3,3) * ones(n,1);
tmp_x(vertex_Cubic_idx_in,1) = anisotropic_mtx(3,3);
i_idx(leng+1:leng+n,1)       = (5*n+1:6*n)';
j_idx(leng+1:leng+n,1)       = (5*n+1:6*n)';
val(leng+1:leng+n,1)         = tmp_x;
leng                         = leng + n;

%
% block (2,3)
%
tmp_x                       = value_air(2,3) * ones(n,1);
tmp_x(face_edge_idx_in.X,1) = anisotropic_mtx(2,3);
i_idx(leng+1:leng+n,1)      = (4*n+1:5*n)';
j_idx(leng+1:leng+n,1)      = (8*n+1:9*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(3,2) * ones(n,1);
tmp_x(vertex_Cubic_idx_in,1) = anisotropic_mtx(3,2);
i_idx(leng+1:leng+n,1)      = (5*n+1:6*n)';
j_idx(leng+1:leng+n,1)      = (7*n+1:8*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

% 
% block (2,4)
%
tmp_x                       = value_air(1,3) * ones(n,1);
tmp_x(face_edge_idx_in.Y,1) = anisotropic_mtx(1,3);
i_idx(leng+1:leng+n,1)      = (3*n+1:4*n)';
j_idx(leng+1:leng+n,1)      = (11*n+1:12*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(3,1) * ones(n,1);
tmp_x(vertex_Cubic_idx_in,1) = anisotropic_mtx(3,1);
i_idx(leng+1:leng+n,1)      = (5*n+1:6*n)';
j_idx(leng+1:leng+n,1)      = (9*n+1:10*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

% ======================================================
%
% block (3,1)
%
tmp_x                       = value_air(1,3) * ones(n,1);
tmp_x(face_edge_idx_in.Z,1) = anisotropic_mtx(1,3);
i_idx(leng+1:leng+n,1)      = (6*n+1:7*n)';
j_idx(leng+1:leng+n,1)      = (2*n+1:3*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(3,1) * ones(n,1);
tmp_x(face_edge_idx_in.X,1) = anisotropic_mtx(3,1);
i_idx(leng+1:leng+n,1)      = (8*n+1:9*n)';
j_idx(leng+1:leng+n,1)      = (1:n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

%
% block (3,2)
%
tmp_x                       = value_air(2,3) * ones(n,1);
tmp_x(vertex_Cubic_idx_in,1) = anisotropic_mtx(2,3);
i_idx(leng+1:leng+n,1)      = (7*n+1:8*n)';
j_idx(leng+1:leng+n,1)      = (5*n+1:6*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(3,2) * ones(n,1);
tmp_x(face_edge_idx_in.X,1) = anisotropic_mtx(3,2);
i_idx(leng+1:leng+n,1)      = (8*n+1:9*n)';
j_idx(leng+1:leng+n,1)      = (4*n+1:5*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

%
% block (3,3)
%
tmp_x                       = value_air(1,1) * ones(n,1);
tmp_x(face_edge_idx_in.Z,1) = anisotropic_mtx(1,1);
i_idx(leng+1:leng+n,1)      = (6*n+1:7*n)';
j_idx(leng+1:leng+n,1)      = (6*n+1:7*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(2,2) * ones(n,1);
tmp_x(vertex_Cubic_idx_in,1) = anisotropic_mtx(2,2);
i_idx(leng+1:leng+n,1)      = (7*n+1:8*n)';
j_idx(leng+1:leng+n,1)      = (7*n+1:8*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(3,3) * ones(n,1);
tmp_x(face_edge_idx_in.X,1) = anisotropic_mtx(3,3);
i_idx(leng+1:leng+n,1)      = (8*n+1:9*n)';
j_idx(leng+1:leng+n,1)      = (8*n+1:9*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

%
% block (3,4)
%

tmp_x                       = value_air(1,2) * ones(n,1);
tmp_x(face_edge_idx_in.Z,1) = anisotropic_mtx(1,2);
i_idx(leng+1:leng+n,1)      = (6*n+1:7*n)';
j_idx(leng+1:leng+n,1)      = (10*n+1:11*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(2,1) * ones(n,1);
tmp_x(vertex_Cubic_idx_in,1) = anisotropic_mtx(2,1);
i_idx(leng+1:leng+n,1)      = (7*n+1:8*n)';
j_idx(leng+1:leng+n,1)      = (9*n+1:10*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

% ======================================================
%
% block (4,1)
%
tmp_x                       = value_air(2,3) * ones(n,1);
tmp_x(face_edge_idx_in.Z,1) = anisotropic_mtx(2,3);
i_idx(leng+1:leng+n,1)      = (10*n+1:11*n)';
j_idx(leng+1:leng+n,1)      = (2*n+1:3*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(3,2) * ones(n,1);
tmp_x(face_edge_idx_in.Y,1) = anisotropic_mtx(3,2);
i_idx(leng+1:leng+n,1)      = (11*n+1:12*n)';
j_idx(leng+1:leng+n,1)      = (n+1:2*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

%
% block (4,2)
%
tmp_x                       = value_air(1,3) * ones(n,1);
tmp_x(vertex_Cubic_idx_in,1) = anisotropic_mtx(1,3);
i_idx(leng+1:leng+n,1)      = (9*n+1:10*n)';
j_idx(leng+1:leng+n,1)      = (5*n+1:6*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(3,1) * ones(n,1);
tmp_x(face_edge_idx_in.Y,1) = anisotropic_mtx(3,1);
i_idx(leng+1:leng+n,1)      = (11*n+1:12*n)';
j_idx(leng+1:leng+n,1)      = (3*n+1:4*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

%
% block (4,3)
%
tmp_x                       = value_air(1,2) * ones(n,1);
tmp_x(vertex_Cubic_idx_in,1) = anisotropic_mtx(1,2);
i_idx(leng+1:leng+n,1)      = (9*n+1:10*n)';
j_idx(leng+1:leng+n,1)      = (7*n+1:8*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(2,1) * ones(n,1);
tmp_x(face_edge_idx_in.Z,1) = anisotropic_mtx(2,1);
i_idx(leng+1:leng+n,1)      = (10*n+1:11*n)';
j_idx(leng+1:leng+n,1)      = (6*n+1:7*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

%
% block (4,4)
%
tmp_x                       = value_air(1,1) * ones(n,1);
tmp_x(vertex_Cubic_idx_in,1) = anisotropic_mtx(1,1);
i_idx(leng+1:leng+n,1)      = (9*n+1:10*n)';
j_idx(leng+1:leng+n,1)      = (9*n+1:10*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(2,2) * ones(n,1);
tmp_x(face_edge_idx_in.Z,1) = anisotropic_mtx(2,2);
i_idx(leng+1:leng+n,1)      = (10*n+1:11*n)';
j_idx(leng+1:leng+n,1)      = (10*n+1:11*n)';
val(leng+1:leng+n,1)        = tmp_x;
leng                        = leng + n;

tmp_x                       = value_air(3,3) * ones(n,1);
tmp_x(face_edge_idx_in.Y,1) = anisotropic_mtx(3,3);
i_idx(leng+1:leng+n,1)      = (11*n+1:12*n)';
j_idx(leng+1:leng+n,1)      = (11*n+1:12*n)';
val(leng+1:leng+n,1)        = tmp_x;
%leng                        = leng + n;

mtx_parameter       = sparse(i_idx, j_idx, val);

[i_idx, j_idx, val] = find(mtx_parameter);
mtx_parameter       =sparse(i_idx, j_idx, val);

end

