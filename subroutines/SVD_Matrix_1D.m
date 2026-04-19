function [Lambda, V, U] = SVD_Matrix_1D(n, m, len_cell, BC_info)
    global U_m V_m F_m ErrValidation flag_numerical_validation 
    
    switch BC_info.str
        case {'Periodic','periodic'}  
            Lambda = spdiags( DistOp.lambda_m.Periodic(n,len_cell),...
                                 0, n, m);
            Lambda = Lambda(2:end,2:end);
            
            V.R = FAME_fnc.R_m_fnc.DFT(n); 
            V.P = FAME_fnc.P_m_fnc.DFT(n); 
            V.D = FAME_fnc.D_m_fnc.DFT(n); 
            V.start_idx = 1;
            U.R = FAME_fnc.R_m_fnc.DFT(n); 
            U.P = FAME_fnc.P_m_fnc.DFT(n); 
            U.D = FAME_fnc.D_m_fnc.DFT(n);
            U.start_idx = 1;
            if flag_numerical_validation
                K_m = FAME_fnc.K_m_PBC(n,len_cell);
                F_m{end+1} = FAME_fnc.F_m_fnc(n);
                U_m{end+1} = U.R*F_m{end}*U.P;
                V_m{end+1} = V.R*F_m{end}*V.P;
                
                x = rand(m,1) + 1i*rand(m,1);
                ErrValidation.FFT_F_m(end+1) = norm(ifft(x)*sqrt(m) - F_m{end}*x);
                ErrValidation.DFT_V_m(end+1) = -1;
                ErrValidation.DFT_U_m(end+1) = -1;
            end
        case {'Quasi-Periodic','QP','quasi-periodic','quasiperiodic'}  
            Lambda = spdiags( DistOp.lambda_QuasiPeriodic(n,BC_info.theta,len_cell),...
                                 0, n, m);
            
            V.type = 'DFT';
            V.R = DistOp.R_sDFT(n);
            V.P = DistOp.P_sDFT(n); 
            V.D = DistOp.D_sDFT(n);
            V.S = DistOp.S_sDFT(n,BC_info.theta);
            V.P_idx = DistOp.Pidx_DFT(n);
            [V.R_idx1, V.R_idx2, V.R_sgn] = DistOp.Ridx_DFT(n);
            V.start_idx = 1;
            U.type = 'DFT';
            U.R = DistOp.R_sDFT(n); 
            U.P = DistOp.P_sDFT(n); 
            U.D = DistOp.D_sDFT(n);
            U.S = DistOp.S_sDFT(n,BC_info.theta);
            U.P_idx = DistOp.Pidx_DFT(n);
            [U.R_idx1,U.R_idx2, U.R_sgn] = DistOp.Ridx_DFT(n);
            U.start_idx = 1;
            if flag_numerical_validation
                K_m = DistOp.K_QPBC(n,BC_info.theta,len_cell);
                F_m{end+1} = DistOp.DFT_mtx(n);
                U_m{end+1} = diag(V.S)*F_m{end};
                V_m{end+1} = diag(V.S)*F_m{end};
                
                x = rand(m,1) + 1i*rand(m,1);
                ErrValidation.FFT_F_m(end+1) = norm(ifft(x)*sqrt(m) - F_m{end}*x);
                ErrValidation.DFT_V_m(end+1) = -1;
                ErrValidation.DFT_U_m(end+1) = -1;
            end
                             
        case {'Dirichlet','dirichlet','dirich','PEC'}
            Lambda = spdiags( DistOp.lambda_Dirichlet(n,len_cell),...
                                -1, n, m);
            V.type = 'DCT';
            V.R = DistOp.R_DCT(n); 
            V.P = DistOp.P_DCT(n); 
            V.D = DistOp.D_DCT(n);
            V.S = DistOp.S_DCT(n);
            V.P_idx = DistOp.Pidx_DCT(n);
            [V.R_idx1, V.R_idx2, V.R_sgn] = DistOp.Ridx_DCT(n);
            V.start_idx = 1;
            U.type = 'DST';
            U.R = DistOp.R_DST(n); 
            U.P = DistOp.P_DST(n); 
            U.D = DistOp.D_DST(n);
            U.S = DistOp.S_DST(n);
            U.P_idx = DistOp.Pidx_DST(n);
            [U.R_idx1, U.R_idx2, U.R_sgn] = DistOp.Ridx_DST(n);
            U.start_idx = 2;
            if flag_numerical_validation
                K_m  = DistOp.K_DBC(n,len_cell);
                U_m{end+1}  = DistOp.DST_mtx(m);
                V_m{end+1}  = DistOp.DCT_mtx(n);
                F_m{end+1}  = DistOp.DFT_mtx(2*n);
                
                ErrValidation.FFT_F_m(end+1) = -1;
                ErrValidation.DFT_V_m(end+1) = norm(V_m{end} - V.R*F_m{end}*V.P*spdiags(V.D,0,n,n),'fro');
                ErrValidation.DFT_U_m(end+1) = norm(U_m{end} - U.R*F_m{end}*U.P*spdiags(U.D,0,m,m),'fro');
            end
    end
    if flag_numerical_validation
        ErrValidation.SVD_K_m(end+1) = norm(K_m - V_m{end}*Lambda*U_m{end}','fro');
    end
end