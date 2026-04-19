%% DistOp
%  A class of functions for discrete operators (matrices or vectors).
%
%% Contribution
%  Author : Jia-Wei Lin (jiaweilin@nycu.edu.tw)
%  Created: 2022/03/16
% 
%  Copyright Jia-Wei Lin
%  https://sites.google.com/g2.nctu.edu.tw/jwlin
classdef DistOp    
    methods (Static)
        val = theta(m,k); % phase used by SVDs
        
        vec = lambda_Dirichlet(m,delta) % singular values of m x (m-1) K_DBC
        vec = lambda_QuasiPeriodic(m,theta,delta) % singular values of m x m K_QPBC
        vec = lambda_Periodic(m,delta) % singular values of m x m K_PBC
        
        mtx = eye_left(m,n) % left-part of an m x m identity matrix
        mtx = eye_right(m,n) % right-part of an m x m identity matrix
        mtx = eye_flip(m) % flip of m x m identity matrix
        
        mtx = K_DBC(m,delta) % m-by-(m-1) forward difference matrix with DBC
        mtx = K_QPBC(m,theta,delta) % m-by-m forward difference matrix with QPBC
        mtx = K_PBC(m,delta) % m-by-m forward difference matrix with PBC
        mtx = DFT_mtx(m) % m-by-m DFT matrix
        mtx = DCT_mtx(m) % m-by-m DCT matrix
        mtx = DST_mtx(m) % m-by-m DST matrix
        
        % Factorizations of DFT, DCT, and DST matrices
        mtx = P_DFT(m)
        mtx = R_DFT(m)
        vec = D_DFT(m)
        
        vec = S_sDFT(m,theta)
        mtx = P_sDFT(m)
        mtx = R_sDFT(m)
        vec = D_sDFT(m)
        vec = Pidx_DFT(m)
        [vec1, vec2, sgn] = Ridx_DFT(m)
        
        vec = S_DCT(m)
        mtx = P_DCT(m)
        mtx = R_DCT(m)
        vec = D_DCT(m)
        vec = Pidx_DCT(m)
        [vec1, vec2, sgn] = Ridx_DCT(m)
        
        vec = S_DST(m)
        mtx = P_DST(m)
        mtx = R_DST(m)
        vec = D_DST(m)
        vec = Pidx_DST(m)
        [vec1, vec2, sgn] = Ridx_DST(m)
        
        % Curl operator for approximate electric field e
        mtx = K(m,delta,BC_info) % forward difference matrix
        mtx = C_32(row,col,delta,BC) % y-partial derivative  for e1 component
        mtx = C_23(row,col,delta,BC) % z-partial derivative  for e1 component
        mtx = C_31(row,col,delta,BC) % x-partial derivative  for e2 component
        mtx = C_13(row,col,delta,BC) % z-partial derivative  for e2 component
        mtx = C_21(row,col,delta,BC) % x-partial derivative  for e3 component
        mtx = C_12(row,col,delta,BC) % y-partial derivative  for e3 component
        mtx = O_x(row,col) % zero matrix
        mtx = O_y(row,col) % zero matrix
        mtx = O_z(row,col) % zero matrix
        mtx = C(row,col,delta,BC) % Curl operator for approximate electric field e
        
        % partial derivative for D and B
        mtx = C_d1(row,col,delta,BC)
        mtx = C_d2(row,col,delta,BC)
        mtx = C_d3(row,col,delta,BC)
        mtx = C_b1(row,col,delta,BC)
        mtx = C_b2(row,col,delta,BC)
        mtx = C_b3(row,col,delta,BC)
        
        mtx = I_R1(row,col)
        mtx = I_R2(row,col)
        mtx = I_R3(row,col)
        mtx = I_S1(row,col)
        mtx = I_S2(row,col)
        mtx = I_S3(row,col)
        
        % Permutation matrices for SVD of C
        mtx = I_Psi_1(row,col)
        mtx = I_Psi_2(row,col)
        mtx = I_Psi_3(row,col)
        mtx = I_Psi(row,col)
        mtx = I_Phi_1(row,col)
        mtx = I_Phi_2(row,col)
        mtx = I_Phi_3(row,col)
        mtx = I_Phi(row,col)
        mtx = I_hat_1(row,col)
        mtx = I_hat_2(row,col)
        mtx = I_hat_3(row,col)
        mtx = I_hat(row,col)
        mtx = I_tilde_1(row,col)
        mtx = I_tilde_2(row,col)
        mtx = I_tilde_3(row,col)
        mtx = I_tilde(row,col)
        
        % SVD of K
%         [mtx_Lambda, mtx_V, mtx_U, val_row, val_col ] = SVD_Matrix_1D(num_cell, len_cell, BC_info)
%         [mtx_Lambda_3D, mtx_R, mtx_P, vec_D, Dim] = SVD_Matrix_3D(Lambda_1D, V, U, Dim)
%         [mtx_Lambda, mtx_Pi, mtx_I] = SVD_Matrix_curl(Lambda_1D, Lambda_3D, Dim)
        
        % matrix-vector production of Psi and Phi
%         vec = MtxVecMul_Psi(x, Dim, type, R, P, D)
%         vec = MtxVecMul_Phi(x, Dim, type, R, P, D)
    end
end