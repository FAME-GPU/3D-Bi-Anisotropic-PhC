function [mtx_Eps_edge, mtx_Eps_face, mtx_Mu_edge, mtx_Mu_face, mtx_Zeta_edge, mtx_Zeta_face, mtx_Xi_edge, mtx_Xi_face] = ...
    mtx_Eps_Mu_Zeta_Xi(n, edge_idx_in, Cubic_idx_in, face_idx_in, vertex_idx_in, mtx_par_eps, mtx_par_mu, mtx_par_zeta, ...
    mtx_par_xi, val_eps_air, val_mu_air, val_zeta_air, val_xi_air)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Pi = sparse(12 * n, 12 * n);
% Pi(1 : n, 1 : n) = eye(n);
% Pi(3 * n + 1 : 4 * n, 3 * n + 1 : 4 * n) = eye(n);
% Pi(6 * n + 1 : 7 * n, 6 * n + 1 : 7 * n) = eye(n);
% Pi(9 * n + 1 : 10 * n, 9 * n + 1 : 10 * n) = eye(n);
% Pi(n + 1 : 2 * n, 4 * n + 1 : 5 * n) = eye(n);
% Pi(4 * n + 1 : 5 * n, n + 1 : 2 * n) = eye(n);
% Pi(7 * n + 1 : 8 * n, 10 * n + 1 : 11 * n) = eye(n);
% Pi(10 * n + 1 : 11 * n, 7 * n + 1 : 8 * n) = eye(n);
% Pi(2 * n + 1 : 3 * n, 8 * n + 1 : 9 * n) = eye(n);
% Pi(8 * n + 1 : 9 * n, 2 * n + 1 : 3 * n) = eye(n);
% Pi(5 * n + 1 : 6 * n, 11 * n + 1 : 12 * n) = eye(n);
% Pi(11 * n + 1 : 12 * n, 5 * n + 1 : 6 * n) = eye(n);
% mtx_par_eps = rand(3);

mtx_Eps_edge = Construct_parameter_mtx(n, edge_idx_in, Cubic_idx_in, mtx_par_eps, val_eps_air);
% Eps_e = Pi' * mtx_Eps_edge * Pi;
% spy(Eps_e)
mtx_Eps_face = Construct_parameter_mtx(n, face_idx_in, vertex_idx_in, mtx_par_eps, val_eps_air);

mtx_Mu_edge = Construct_parameter_mtx(n, edge_idx_in, Cubic_idx_in, mtx_par_mu, val_mu_air);
mtx_Mu_face = Construct_parameter_mtx(n, face_idx_in, vertex_idx_in, mtx_par_mu, val_mu_air);

mtx_Zeta_edge = Construct_parameter_mtx(n, edge_idx_in, Cubic_idx_in, mtx_par_zeta, val_zeta_air);
mtx_Zeta_face = Construct_parameter_mtx(n, face_idx_in, vertex_idx_in, mtx_par_zeta, val_zeta_air);

mtx_Xi_edge = Construct_parameter_mtx(n, edge_idx_in, Cubic_idx_in, mtx_par_xi, val_xi_air);
mtx_Xi_face = Construct_parameter_mtx(n, face_idx_in, vertex_idx_in, mtx_par_xi, val_xi_air);

end

