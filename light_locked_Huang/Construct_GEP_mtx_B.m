function [GEP_mtx_B] = Construct_GEP_mtx_B(n, mtx_Eps_edge, mtx_Eps_face, mtx_Mu_edge, mtx_Mu_face, mtx_Zeta_edge, ...
    mtx_Zeta_face, mtx_Xi_edge, mtx_Xi_face)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[i_idx_Eps_edge, j_idx_Eps_edge, val_idx_Eps_edge ] = find(mtx_Eps_edge);
[i_idx_Eps_face, j_idx_Eps_face, val_idx_Eps_face ] = find(mtx_Eps_face);

[i_idx_Mu_edge, j_idx_Mu_edge, val_idx_Mu_edge ] = find(mtx_Mu_edge);
[i_idx_Mu_face, j_idx_Mu_face, val_idx_Mu_face ] = find(mtx_Mu_face);

[i_idx_Zeta_edge, j_idx_Zeta_edge, val_idx_Zeta_edge ] = find(mtx_Zeta_edge);
[i_idx_Zeta_face, j_idx_Zeta_face, val_idx_Zeta_face ] = find(mtx_Zeta_face);

[i_idx_Xi_edge, j_idx_Xi_edge, val_idx_Xi_edge ] = find(mtx_Xi_edge);
[i_idx_Xi_face, j_idx_Xi_face, val_idx_Xi_face ] = find(mtx_Xi_face);

vec_leng = [ length(i_idx_Eps_edge); length(i_idx_Eps_face); length(i_idx_Mu_edge); length(i_idx_Mu_face); ...
    length(i_idx_Zeta_edge); length(i_idx_Zeta_face); length(i_idx_Xi_edge); length(i_idx_Xi_face) ];

mm = sum(vec_leng);

i_idx = zeros(mm,1);
j_idx = zeros(mm,1);
val   = zeros(mm,1);
%
i_idx(1:vec_leng(6,1),1) = i_idx_Zeta_face;
j_idx(1:vec_leng(6,1),1) = j_idx_Zeta_face + 12 * n;
val(1:vec_leng(6,1),1)   = val_idx_Zeta_face;

leng = vec_leng(6,1);
i_idx(leng+1:leng+vec_leng(5,1),1) = i_idx_Zeta_edge + 12 * n;
j_idx(leng+1:leng+vec_leng(5,1),1) = j_idx_Zeta_edge;
val(leng+1:leng+vec_leng(5,1),1)   = val_idx_Zeta_edge;
%
leng = leng + vec_leng(5,1);
i_idx(leng+1:leng+vec_leng(4,1),1) = i_idx_Mu_face;
j_idx(leng+1:leng+vec_leng(4,1),1) = j_idx_Mu_face + 24 * n;
val(leng+1:leng+vec_leng(4,1),1)   = val_idx_Mu_face;

leng = leng + vec_leng(4,1);
i_idx(leng+1:leng+vec_leng(3,1),1) = i_idx_Mu_edge + 12 * n;
j_idx(leng+1:leng+vec_leng(3,1),1) = j_idx_Mu_edge + 36 * n;
val(leng+1:leng+vec_leng(3,1),1)   = val_idx_Mu_edge;
%
leng = leng + vec_leng(3,1);
i_idx(leng+1:leng+vec_leng(1,1),1) = i_idx_Eps_edge + 24 * n;
j_idx(leng+1:leng+vec_leng(1,1),1) = j_idx_Eps_edge;
val(leng+1:leng+vec_leng(1,1),1)   = -val_idx_Eps_edge;

leng = leng + vec_leng(1,1);
i_idx(leng+1:leng+vec_leng(2,1),1) = i_idx_Eps_face + 36 * n;
j_idx(leng+1:leng+vec_leng(2,1),1) = j_idx_Eps_face + 12 * n;
val(leng+1:leng+vec_leng(2,1),1)   = -val_idx_Eps_face;
%
leng = leng + vec_leng(2,1);
i_idx(leng+1:leng+vec_leng(7,1),1) = i_idx_Xi_edge + 24 * n;
j_idx(leng+1:leng+vec_leng(7,1),1) = j_idx_Xi_edge + 36 * n;
val(leng+1:leng+vec_leng(7,1),1)   = -val_idx_Xi_edge;

leng = leng + vec_leng(7,1);
i_idx(leng+1:leng+vec_leng(8,1),1) = i_idx_Xi_face + 36 * n;
j_idx(leng+1:leng+vec_leng(8,1),1) = j_idx_Xi_face + 24 * n;
val(leng+1:leng+vec_leng(8,1),1)   = -val_idx_Xi_face;

GEP_mtx_B = sparse(i_idx, j_idx, val);

end

