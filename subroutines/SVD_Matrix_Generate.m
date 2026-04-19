function [ mtx ] = SVD_Matrix_Generate(opt)
    global V_m U_m F_m 
    global flag_numerical_validation ErrValidation C
    flag_numerical_validation = opt.comp.flag_numerical_validation;
    ErrValidation = {};
    if flag_numerical_validation
        V_m = {}; U_m = {}; F_m = {};
        ErrValidation.SVD_K_m = [];   
        ErrValidation.DFT_V_m = [];   
        ErrValidation.DFT_U_m = [];   
        ErrValidation.FFT_F_m = []; 
    end
    
    hx = opt.comp.discrete.mesh_len(1);
    hy = opt.comp.discrete.mesh_len(2);
    hz = opt.comp.discrete.mesh_len(3);
    
    Dim = opt.comp.Dim;
    n1 = Dim.n1; m1 = Dim.m1;
    n2 = Dim.n2; m2 = Dim.m2;
    n3 = Dim.n3; m3 = Dim.m3;    
    %% Construct SVD matrices
    % 1D singular value matrices
    [Lambda_1D.x, mtx.V.x, mtx.U.x] = SVD_Matrix_1D(n1, m1, hx, opt.phys.BC.x);
    [Lambda_1D.y, mtx.V.y, mtx.U.y] = SVD_Matrix_1D(n2, m2, hy, opt.phys.BC.y);
    [Lambda_1D.z, mtx.V.z, mtx.U.z] = SVD_Matrix_1D(n3, m3, hz, opt.phys.BC.z);
    % 3D singular value matrices of discrete partial derivatives C_1, C_2, and C_3
    [Lambda_3D, mtx.R, mtx.R_idx1, mtx.R_idx2, mtx.P, mtx.P_idx, mtx.D, mtx.S, mtx.Type] = SVD_Matrix_3D(Lambda_1D, mtx.V, mtx.U, Dim);
    % 3D singular value matrices of discrete curl C
    [mtx.Lambda_curl,mtx.Pi,mtx.I] = SVD_Matrix_curl(Lambda_1D, Lambda_3D, Dim);

    
    
    
    
    if opt.comp.flag_numerical_validation
        delta.x = hx; delta.y = hy; delta.z = hz;    
        row.x = n1; row.y = n2; row.z = n3;
        col.x = m1; col.y = m2; col.z = m3;
        BC = opt.phys.BC;

        C_1.y = DistOp.C_32(row,col,delta,BC);
        C_1.z = DistOp.C_23(row,col,delta,BC);
        C_2.x = DistOp.C_31(row,col,delta,BC);
        C_2.z = DistOp.C_13(row,col,delta,BC);
        C_3.x = DistOp.C_21(row,col,delta,BC);
        C_3.y = DistOp.C_12(row,col,delta,BC);
        C     = DistOp.C(row,col,delta,BC);
        
        Psi_1 = kron( V_m{3}, kron( V_m{2}, U_m{1} ) );
        Psi_2 = kron( V_m{3}, kron( U_m{2}, V_m{1} ) );
        Psi_3 = kron( U_m{3}, kron( V_m{2}, V_m{1} ) );
        Phi_1 = kron( U_m{3}, kron( U_m{2}, V_m{1} ) );
        Phi_2 = kron( U_m{3}, kron( V_m{2}, U_m{1} ) );
        Phi_3 = kron( V_m{3}, kron( U_m{2}, U_m{1} ) );
        
        Psi   = blkdiag(Psi_1,Psi_2,Psi_3);
        Phi   = blkdiag(Phi_1,Phi_2,Phi_3);
        ErrValidation.SVD_C_2x = norm(C_2.x - Psi_3*Lambda_3D.C31*Phi_2','fro');
        ErrValidation.SVD_C_3x = norm(C_3.x - Psi_2*Lambda_3D.C21*Phi_3','fro');
        ErrValidation.SVD_C_1y = norm(C_1.y - Psi_3*Lambda_3D.C32*Phi_1','fro');
        ErrValidation.SVD_C_3y = norm(C_3.y - Psi_1*Lambda_3D.C12*Phi_3','fro');            
        ErrValidation.SVD_C_1z = norm(C_1.z - Psi_2*Lambda_3D.C23*Phi_1','fro');
        ErrValidation.SVD_C_2z = norm(C_2.z - Psi_1*Lambda_3D.C13*Phi_2','fro');

        ErrValidation.Compact_Lambda_1 = norm(mtx.Lambda_curl.r - mtx.I.Psi'*mtx.Lambda_curl.orig*mtx.I.Phi,'fro');
        ErrValidation.Compact_Lambda_2 = norm(mtx.Lambda_curl.tilde_q - mtx.I.tilde'*mtx.Lambda_curl.orig'*mtx.Lambda_curl.orig*mtx.I.tilde,'fro');
        ErrValidation.Compact_Lambda_3 = norm(mtx.Lambda_curl.orig'*mtx.I.hat,'fro');
        
        Psi_r   = Psi * mtx.I.Psi;
        Phi_r   = Phi * mtx.I.Phi;
        ErrValidation.Compact_C = norm(C*Phi_r - Psi_r*mtx.Lambda_curl.r,'fro');
        ErrValidation.SVD_Lambda_r = norm(mtx.Lambda_curl.r - mtx.Pi.Pq*mtx.Lambda_curl.Sigma.q*mtx.Pi.Qq','fro');
        
        Qq = Phi_r*mtx.Pi.Qq; Qqr = Phi_r*mtx.Pi.Qqr;
        Pq = Psi_r*mtx.Pi.Pq; Pqr = Psi_r*mtx.Pi.Pqr;
        Qtilde = Phi*mtx.I.tilde;
        Ptilde = Psi*mtx.Lambda_curl.orig*mtx.I.tilde*mtx.Lambda_curl.Sigma.tilde_q_inv;
        Ptilde_0 = Psi*mtx.Pi.tilde*mtx.I.tilde*mtx.Lambda_curl.Sigma.tilde_q_inv;
        ErrValidation.PSVD_C_q = norm(C*Qq - Pq*mtx.Lambda_curl.Sigma.q,'fro');
        ErrValidation.PSVD_C_tilde_1 = norm(C*Qtilde - Ptilde*mtx.Lambda_curl.Sigma.tilde_q,'fro');
        ErrValidation.PSVD_C_tilde_2 = norm(C'*Ptilde_0,'fro');
        
        ErrValidation.equalities_1 = norm(mtx.I.Psi'*mtx.Lambda_curl.orig*mtx.I.tilde,'fro');
        ErrValidation.equalities_2 = norm(mtx.I.Psi'*mtx.Pi.tilde*mtx.I.tilde,'fro');
        ErrValidation.equalities_3 = norm(mtx.I.Psi'*mtx.I.hat,'fro');
        ErrValidation.equalities_4 = norm(mtx.Lambda_curl.orig'*mtx.Pi.tilde*mtx.I.tilde,'fro');
        ErrValidation.equalities_5 = norm(mtx.I.tilde'*mtx.Lambda_curl.orig'*mtx.Pi.tilde,'fro');
        ErrValidation.equalities_6 = norm(mtx.Lambda_curl.orig'*mtx.I.hat,'fro');
        ErrValidation.equalities_7 = norm(mtx.Pi.tilde'*mtx.I.hat,'fro');
        
        Phat_0 = Psi*mtx.I.hat;
        ErrValidation.PSVD_C_hat = norm(C'*Phat_0,'fro');
        
        SVD_P = [Phat_0, Ptilde_0, Pq, Ptilde];
        SVD_Q = [Qq, Qtilde];
        SVD_Pr = [ Pqr, Ptilde];
        SVD_Qr = [ Qqr, Qtilde];
        ErrValidation.semiorig_P  = norm( SVD_P'*SVD_P - speye(size(SVD_P,2)), 'fro' );
        ErrValidation.semiorig_Q  = norm( SVD_Q'*SVD_Q - speye(size(SVD_Q,2)), 'fro' );
        ErrValidation.semiorig_Pr = norm( SVD_Pr'*SVD_Pr - speye(size(SVD_Pr,2)), 'fro' );
        ErrValidation.semiorig_Qr = norm( SVD_Qr'*SVD_Qr - speye(size(SVD_Qr,2)), 'fro' );
    end
    
end




