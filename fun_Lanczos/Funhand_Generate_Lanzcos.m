function fnc = Funhand_Generate_Lanzcos( Dim, mtx_A, mtx_B, linearsolver, NFSEP_type  )
global flag_numerical_validation C ErrValidation flag_gpu flag_FFT
global time_tmp time_pcg
global time_Psis_tmp time_Psi_tmp time_Phis_tmp time_Phi_tmp
global time_FFT
    
    time_pcg.invtildeM = EventTime.pcg('M');
    time_pcg.invtildeN = EventTime.pcg('N');


    % Matrix vector mutiplication of DFT, DCT and DST corresponding
    % to DFT_mtx, DCT_mtx, and DST_mtx
    fnc.DFT_x = @(x) sqrt(size(x,1))\fft(x, [], 1); fnc.IDFT_x = @(x) sqrt(size(x,1))*ifft(x, [], 1); 
    fnc.DFT_y = @(x) sqrt(size(x,2))\fft(x, [], 2); fnc.IDFT_y = @(x) sqrt(size(x,2))*ifft(x, [], 2); 
    fnc.DFT_z = @(x) sqrt(size(x,3))\fft(x, [], 3); fnc.IDFT_z = @(x) sqrt(size(x,3))*ifft(x, [], 3); 
    
    fnc.DCT_x = @(x) dct(x, [], 1, 'type', 2); fnc.IDCT_x = @(x) idct(x, [], 1, 'type', 2); 
    fnc.DCT_y = @(x) dct(x, [], 2, 'type', 2); fnc.IDCT_y = @(x) idct(x, [], 2, 'type', 2); 
    fnc.DCT_z = @(x) dct(x, [], 3, 'type', 2); fnc.IDCT_z = @(x) idct(x, [], 3, 'type', 2); 
    
    if flag_gpu
        fnc.DST_x = @(x) gDST(real(x), [], 1 ) + 1i*gDST(imag(x), [], 1 ); fnc.IDST_x = @(x) gDST(real(x), [], 1 ) + 1i*gDST(imag(x), [], 1 ); 
        fnc.DST_y = @(x) gDST(real(x), [], 2 ) + 1i*gDST(imag(x), [], 2 ); fnc.IDST_y = @(x) gDST(real(x), [], 2 ) + 1i*gDST(imag(x), [], 2 );
        fnc.DST_z = @(x) gDST(real(x), [], 3 ) + 1i*gDST(imag(x), [], 3 ); fnc.IDST_z = @(x) gDST(real(x), [], 3 ) + 1i*gDST(imag(x), [], 3 );
    else
        fnc.DST_x = @(x) DST(real(x), [], 1 ) + 1i*DST(imag(x), [], 1 ); fnc.IDST_x = @(x) DST(real(x), [], 1 ) + 1i*DST(imag(x), [], 1 ); 
        fnc.DST_y = @(x) DST(real(x), [], 2 ) + 1i*DST(imag(x), [], 2 ); fnc.IDST_y = @(x) DST(real(x), [], 2 ) + 1i*DST(imag(x), [], 2 );
        fnc.DST_z = @(x) DST(real(x), [], 3 ) + 1i*DST(imag(x), [], 3 ); fnc.IDST_z = @(x) DST(real(x), [], 3 ) + 1i*DST(imag(x), [], 3 );
    end
    
%     x = (1:60)';
%     x = reshape(x, 3,4,5);
%     y = fnc.DST_x(x)
    %% Singular vectors of C_pec
    % Construct DFT-based matrix-vector multiplication of Psi and Phi
    if flag_gpu
        
        fnc.Psis = @(x) gMtxVecMul_FFT(x, Dim.Psi,   'transpose', mtx_A.S.Psi, mtx_A.D.Psi, mtx_A.P_idx.Psi_1, mtx_A.P_idx.Psi_2, mtx_A.P_idx.Psi_3, mtx_A.Type.Psi_1, mtx_A.Type.Psi_2, mtx_A.Type.Psi_3);
        fnc.Phis = @(x) gMtxVecMul_FFT(x, Dim.Phi,   'transpose', mtx_A.S.Phi, mtx_A.D.Phi, mtx_A.P_idx.Phi_1, mtx_A.P_idx.Phi_2, mtx_A.P_idx.Phi_3, mtx_A.Type.Phi_1, mtx_A.Type.Phi_2, mtx_A.Type.Phi_3);
        fnc.Psi  = @(x) gMtxVecMul_FFT(x, Dim.Psi, 'notranspose', mtx_A.S.Psi, mtx_A.D.Psi, mtx_A.P_idx.Psi_1, mtx_A.P_idx.Psi_2, mtx_A.P_idx.Psi_3, mtx_A.Type.Psi_1, mtx_A.Type.Psi_2, mtx_A.Type.Psi_3);
        fnc.Phi  = @(x) gMtxVecMul_FFT(x, Dim.Phi, 'notranspose', mtx_A.S.Phi, mtx_A.D.Phi, mtx_A.P_idx.Phi_1, mtx_A.P_idx.Phi_2, mtx_A.P_idx.Phi_3, mtx_A.Type.Phi_1, mtx_A.Type.Phi_2, mtx_A.Type.Phi_3);     
        
    else
        if flag_FFT
            fnc.Psis = @(x) fnc_Psis(x);
            fnc.Phis = @(x) fnc_Phis(x); %MtxVecMul_FFT(x, Dim.Phi,   'transpose', mtx_A.S.Phi, mtx_A.D.Phi, mtx_A.P_idx.Phi_1, mtx_A.P_idx.Phi_2, mtx_A.P_idx.Phi_3, mtx_A.Type.Phi_1, mtx_A.Type.Phi_2, mtx_A.Type.Phi_3);
            fnc.Psi  = @(x) fnc_Psi(x);  %MtxVecMul_FFT(x, Dim.Psi, 'notranspose', mtx_A.S.Psi, mtx_A.D.Psi, mtx_A.P_idx.Psi_1, mtx_A.P_idx.Psi_2, mtx_A.P_idx.Psi_3, mtx_A.Type.Psi_1, mtx_A.Type.Psi_2, mtx_A.Type.Psi_3);
            fnc.Phi  = @(x) fnc_Phi(x);  %MtxVecMul_FFT(x, Dim.Phi, 'notranspose', mtx_A.S.Phi, mtx_A.D.Phi, mtx_A.P_idx.Phi_1, mtx_A.P_idx.Phi_2, mtx_A.P_idx.Phi_3, mtx_A.Type.Phi_1, mtx_A.Type.Phi_2, mtx_A.Type.Phi_3);     
        else
            % Construct DFT/DCT/DST based matrix-vector multiplication of Psi and Phi
            fnc.Psis = @(x) MtxVecMul_DFCST(x, Dim.H1, Dim.H2, Dim.H3, mtx_A.Type.Psi_1, mtx_A.Type.Psi_2, mtx_A.Type.Psi_3,   'transpose', fnc, mtx_A.S.Psi);
            fnc.Psi  = @(x) MtxVecMul_DFCST(x, Dim.H1, Dim.H2, Dim.H3, mtx_A.Type.Psi_1, mtx_A.Type.Psi_2, mtx_A.Type.Psi_3, 'notranspose', fnc, mtx_A.S.Psi);
            fnc.Phis = @(x) MtxVecMul_DFCST(x, Dim.E1, Dim.E2, Dim.E3, mtx_A.Type.Phi_1, mtx_A.Type.Phi_2, mtx_A.Type.Phi_3,   'transpose', fnc, mtx_A.S.Phi);
            fnc.Phi  = @(x) MtxVecMul_DFCST(x, Dim.E1, Dim.E2, Dim.E3, mtx_A.Type.Phi_1, mtx_A.Type.Phi_2, mtx_A.Type.Phi_3, 'notranspose', fnc, mtx_A.S.Phi);
        end
    end
    function y = fnc_Psis(x)
        time_tmp = time_Psis_tmp;
        y = MtxVecMul_FFT(x, Dim.Psi,   'transpose', mtx_A.S.Psi, mtx_A.D.Psi, mtx_A.P_idx.Psi_1, mtx_A.P_idx.Psi_2, mtx_A.P_idx.Psi_3, mtx_A.Type.Psi_1, mtx_A.Type.Psi_2, mtx_A.Type.Psi_3);
        time_Psis_tmp = time_tmp;
    end
    function y = fnc_Psi(x)
        time_tmp = time_Psi_tmp;
        y = MtxVecMul_FFT(x, Dim.Psi, 'notranspose', mtx_A.S.Psi, mtx_A.D.Psi, mtx_A.P_idx.Psi_1, mtx_A.P_idx.Psi_2, mtx_A.P_idx.Psi_3, mtx_A.Type.Psi_1, mtx_A.Type.Psi_2, mtx_A.Type.Psi_3);
        time_Psi_tmp = time_tmp;
    end
    function y = fnc_Phis(x)
        time_tmp = time_Phis_tmp;
        y = MtxVecMul_FFT(x, Dim.Phi,   'transpose', mtx_A.S.Phi, mtx_A.D.Phi, mtx_A.P_idx.Phi_1, mtx_A.P_idx.Phi_2, mtx_A.P_idx.Phi_3, mtx_A.Type.Phi_1, mtx_A.Type.Phi_2, mtx_A.Type.Phi_3);
        time_Phis_tmp = time_tmp;
    end
    function y = fnc_Phi(x)
        time_tmp = time_Phi_tmp;
        y = MtxVecMul_FFT(x, Dim.Phi, 'notranspose', mtx_A.S.Phi, mtx_A.D.Phi, mtx_A.P_idx.Phi_1, mtx_A.P_idx.Phi_2, mtx_A.P_idx.Phi_3, mtx_A.Type.Phi_1, mtx_A.Type.Phi_2, mtx_A.Type.Phi_3);     
        time_Phi_tmp = time_tmp;
    end
    
    %% Singular Value Decomposition of C
    fnc.Pr  = @(x) fnc.Psi( mtx_A.Pi.Pr*x );
    fnc.Qr  = @(x) fnc.Phi( mtx_A.Pi.Qr*x );
    fnc.Prs = @(x) mtx_A.Pi.Pr'*fnc.Psis( x );
    fnc.Qrs = @(x) mtx_A.Pi.Qr'*fnc.Phis( x );
        
    %% Matrics for NFSEP
    fnc.tildeN     = @(x) fnc.Qrs(mtx_B.N_int*fnc.Qr(x));
    fnc.tildeM     = @(x) fnc.Prs(mtx_B.M_int*fnc.Pr(x));

    fnc.invN       = @(x) pcg1(mtx_B.N_int, x, linearsolver.tol, linearsolver.itmax, [], [], linearsolver.x0_fnc(length(x)));
    fnc.invM       = @(x) pcg1(mtx_B.M_int, x, linearsolver.tol, linearsolver.itmax, [], [], linearsolver.x0_fnc(length(x)));
    fnc.invtildeN  = @(x) fnc_invtildeN(x);
    fnc.invtildeM  = @(x) fnc_invtildeM(x);
    
    fnc.invN       = @(x) pcg1(mtx_B.N_int, x, linearsolver.tol, linearsolver.itmax, [], [], linearsolver.x0_fnc(length(x)));
    fnc.invM       = @(x) pcg1(mtx_B.M_int, x, linearsolver.tol, linearsolver.itmax, [], [], linearsolver.x0_fnc(length(x)));
    %fnc.phi        = @(x) mtx_B.N_int*x - mtx_B.J_int * (fnc.invM(mtx_B.L_int * x));
    fnc.invphi     = @(x) pcg1(fnc.phi, x, linearsolver.tol, linearsolver.itmax, [], [], linearsolver.x0_fnc(length(x)));
    fnc.tilde_Mu   = @(x) fnc.Prs(fnc.invM( fnc.Pr(x)) );
    % fnc.tilde_Mu1   = @(x) fnc.invM( fnc.Pr(x)) ;
    fnc.tilde_phi  = @(x) fnc.Qrs(fnc.invN( fnc.Qr(x)) );       
    fnc.invtilde_Mu   = @(x) pcg1(fnc.tilde_Mu, x, linearsolver.tol, linearsolver.itmax, [], [], linearsolver.x0_fnc(length(x)));
    fnc.invtilde_phi  = @(x) pcg1(fnc.tilde_phi, x, linearsolver.tol, linearsolver.itmax, [], [], linearsolver.x0_fnc(length(x)));
    % load('b.mat');
    % x = rand(128,1);
    % xx1 = fnc.tilde_phi(b);
    % xx2 = fnc.invM(xx1);
    % xx3 = fnc.invM(xx1);
    % xx4 = fnc.invtilde_phi(x);
   % v0 = ones(length(mtx_A.Lambda_curl.Sigma.r),1);
   % v1 = fnc.invphi(v0);
    function y = fnc_invtildeN(x)
        time_Phis_tmp = EventTime.Phi('t');
        time_Phi_tmp  = EventTime.Phi('nt');
        pcg_tic = tic;    %time_pcg.invtildeN.iter(end+1)
        [ y, ~, ~, time_pcg.invtildeN.iter(end+1) ] = pcg( fnc.tildeN, x, linearsolver.tol, linearsolver.itmax, [], [], linearsolver.x0_fnc(length(x)));    
        time_pcg.invtildeN.time(end+1) = toc(pcg_tic);
        time_pcg.invtildeN.Phis(end+1) = time_Phis_tmp;
        time_pcg.invtildeN.Phi(end+1)  = time_Phi_tmp;
        fprintf('|');
    end
    function y = fnc_invtildeM(x)
        time_Psis_tmp = EventTime.Psi('t');
        time_Psi_tmp  = EventTime.Psi('nt');
        pcg_tic = tic;    
        [ y, ~, ~, time_pcg.invtildeM.iter(end+1) ] = pcg( fnc.tildeM, x, linearsolver.tol, linearsolver.itmax, [], [], linearsolver.x0_fnc(length(x)));    
        time_pcg.invtildeM.time(end+1) = toc(pcg_tic);
        time_pcg.invtildeM.Psis(end+1) = time_Psis_tmp;
        time_pcg.invtildeM.Psi(end+1)  = time_Psi_tmp;
        fprintf('|');
    end

    switch NFSEP_type
        case 0
            fnc.Ar    = @(x) mtx_A.Lambda_curl.Sigma.r.*fnc.tildeN(mtx_A.Lambda_curl.Sigma.r.*fnc.tildeM(x));
            fnc.invAr = @(x) fnc.invtildeM(mtx_A.Lambda_curl.Sigma.r.\fnc.invtildeN(mtx_A.Lambda_curl.Sigma.r.\x));
        case 1
%             fnc.Ar    = @(x) mtx_A.Lambda_curl.Sigma.r.*fnc.tildeN(mtx_A.Lambda_curl.Sigma.r.*fnc.tildeM(x));
            sqrt_Sigma_r = sqrt(mtx_A.Lambda_curl.Sigma.r);
            fnc.invAr = @(x) [ -sqrt_Sigma_r.\fnc.invtildeN(sqrt_Sigma_r.\x(      1:end/2))  ;
                                sqrt_Sigma_r.\fnc.invtildeM(sqrt_Sigma_r.\x(end/2+1:end  )) ];
%             fnc.invAr = @(x) funhand_test(x,sqrt_Sigma_r,fnc.invtildeM,fnc.tildeN);
    end
    
   % fnc.Ar_shift  = @(x,sigma) fnc.Ar(x) - sigma*x;
    % fnc.invAr_shift = @(x,sigma) mtx.Lambda_curl.Sigma.r.\fnc.invBr_eps_shift(mtx.Lambda_curl.Sigma.r.\x,sigma);
    %     fnc.P   = @(x) fnc.Psi( mtx.Pi.P*x );
    %     fnc.Ps  = @(x) mtx.Pi.P'*(fnc.Phss( x ));
    %     fnc.Q   = @(x) fnc.Psi( mtx.Pi.P*x );
    %     fnc.Qs  = @(x) mtx.Pi.P'*(fnc.Psis( x ));
    %     fnc.Q0  = @(x) fnc.Phi( mtx.Pi.Q0*x );
    %     fnc.Q0s = @(x) mtx.Pi.Q0'*fnc.Phis( x );
    %     fnc.QBQ_shift = @(x,sigma) fnc.QBQ(x) - sigma*(mtx.Lambda_curl.Sigma.r.\(mtx.Lambda_curl.Sigma.r.\x));
    %     fnc.invB_eps_r_shift = @(x,sigma) pcg(@(x)fnc.QBQ_shift(x,sigma), x, linearsolver.tol, linearsolver.itmax, [], [], x0);
    %     fnc.invB_eps_r_shift = @(x,sigma) bicgstabl(@(x)fnc.QBQ_shift(x,sigma), x, linearsolver.tol, linearsolver.itmax, [], [], x0);
    %     fnc.Ar_shift  = @(x,sigma) fnc.Ar(x) - sigma*x;
    %     fnc.invAr_shift = @(x,sigma) mtx.Lambda_curl.Sigma.r.\fnc.invB_eps_r_shift(mtx.Lambda_curl.Sigma.r.\x,sigma);
    %% Test area
    if flag_numerical_validation 
        x = rand(Dim.Nh,1) + 1i*rand(Dim.Nh,1);
        norm(fnc.Psis(x) - fnc.Psis_DFCST(x))
        norm(fnc.Psi(x)  - fnc.Psi_DFCST(x) )
        x = rand(Dim.Ne,1) + 1i*rand(Dim.Ne,1);
        norm(fnc.Phis(x) - fnc.Phis_DFCST(x))
        norm(fnc.Phi(x)  - fnc.Phi_DFCST(x) )
        
        fnc.Qq  = @(x) fnc.Phi( mtx_A.I.Phi*mtx_A.Pi.Qq*x );
        fnc.Pq  = @(x) fnc.Psi( mtx_A.I.Psi*mtx_A.Pi.Pq*x );
        fnc.tildeQ   = @(x) fnc.Phi( mtx_A.I.tilde*x );
        fnc.tildeP   = @(x) fnc.Psi( mtx_A.Lambda_curl.orig*mtx_A.I.tilde*mtx_A.Lambda_curl.Sigma.tilde_q_inv*x );
        fnc.tildeP_0 = @(x) fnc.Psi( mtx_A.Pi.tilde*mtx_A.I.tilde*mtx_A.Lambda_curl.Sigma.tilde_q_inv*x );
        fnc.hatP_0   = @(x) fnc.Psi( mtx_A.I.hat*x );
        fnc.Q  = @(x) fnc.Phi( mtx_A.Pi.Q*x );
        fnc.Qs = @(x) mtx_A.Pi.Q'*(fnc.Phis( x ));
        fnc.P  = @(x) fnc.Psi( mtx_A.Pi.P*x );
        fnc.Ps = @(x) mtx_A.Pi.P'*(fnc.Psis( x ));
        %% Validation of Corollary 1. factorization of C with Psi and Phi
        x = rand(Dim.Ne,1) + 1i*rand(Dim.Ne,1); ErrValidation.DFT_C  = norm(C *x - fnc.Psi(mtx_A.Lambda_curl.orig *fnc.Phis(x)));
        x = rand(Dim.Nh,1) + 1i*rand(Dim.Nh,1); ErrValidation.DFT_Cs = norm(C'*x - fnc.Phi(mtx_A.Lambda_curl.orig'*fnc.Psis(x)));

        %% Validation of Lemma 1. factoriztations of Lambda
        ErrValidation.lemma1_1 = norm(mtx_A.Lambda_curl.orig*mtx_A.I.Phi - mtx_A.I.Psi*mtx_A.Lambda_curl.r, 'fro' );
        ErrValidation.lemma1_2 = norm(mtx_A.I.tilde'*mtx_A.Lambda_curl.orig'*mtx_A.Lambda_curl.orig*mtx_A.I.tilde - mtx_A.Lambda_curl.tilde_q,'fro');
        ErrValidation.lemma1_3 = norm(mtx_A.Lambda_curl.orig'*mtx_A.I.hat,'fro');
        %% Validation of Lemma 2. compact factorization of C and C^*
        x = rand(3*Dim.Nq,1) + 1i*rand(3*Dim.Nq,1);
        ErrValidation.cDFT_C  = norm(C *fnc.Phi(mtx_A.I.Phi*x) - fnc.Psi(mtx_A.I.Psi*mtx_A.Lambda_curl.r * x ));
        ErrValidation.cDFT_Cs = norm(C'*fnc.Psi(mtx_A.I.Psi*x) - fnc.Phi(mtx_A.I.Phi*mtx_A.Lambda_curl.r'* x ));
        %% Validation of Lemma 3. SVD of Lambda_r
        ErrValidation.SVD_Lambda_r  = norm(mtx_A.Lambda_curl.r - mtx_A.Pi.Pq*mtx_A.Lambda_curl.Sigma.q*mtx_A.Pi.Qq','fro');
        %% Validation of Corollary 2. partial SVD of C with Q_q and P_q
        x = rand(3*Dim.Nq,1) + 1i*rand(3*Dim.Nq,1);
        ErrValidation.pSVD1_C  = norm(C *fnc.Qq(x) - fnc.Pq(mtx_A.Lambda_curl.Sigma.q*x));
        ErrValidation.pSVD1_Cs = norm(C'*fnc.Pq(x) - fnc.Qq(mtx_A.Lambda_curl.Sigma.q*x));
        %% Validation of Corollary 3. partial SVD of C with tildeQ, tildeP, and tildeP_0           
        x = rand(Dim.Ntil,1) + 1i*rand(Dim.Ntil,1);
        ErrValidation.pSVD2_C(1) = norm(C *fnc.tildeQ(x) - fnc.tildeP(mtx_A.Lambda_curl.Sigma.tilde_q*x));
        ErrValidation.pSVD2_C(2) = norm(C'*fnc.tildeP(x) - fnc.tildeQ(mtx_A.Lambda_curl.Sigma.tilde_q*x));
        ErrValidation.pSVD2_C(3) = norm(C'*fnc.tildeP_0(x));  
        %% Validation of Corollary 4. partial SVD of C with hatP_0         
        x = rand(Dim.Nhat,1) + 1i*rand(Dim.Nhat,1);
        ErrValidation.pSVD3_C  = norm(C'*fnc.hatP_0(x));  
        %% Validation of Lemma 4. null subspaces of Lambda, tildePi, and hatI
        ErrValidation.lemma4(1) = norm(mtx_A.I.Psi'*mtx_A.Lambda_curl.orig*mtx_A.I.tilde,'fro');
        ErrValidation.lemma4(2) = norm(mtx_A.I.Psi'*mtx_A.Pi.tilde*mtx_A.I.tilde,'fro');
        ErrValidation.lemma4(3) = norm(mtx_A.I.Psi'*mtx_A.I.hat,'fro');
        ErrValidation.lemma4(4) = norm(mtx_A.Lambda_curl.orig'*mtx_A.Pi.tilde*mtx_A.I.tilde,'fro');
        ErrValidation.lemma4(5) = norm(mtx_A.I.tilde'*mtx_A.Lambda_curl.orig'*mtx_A.Pi.tilde,'fro');
        ErrValidation.lemma4(6) = norm(mtx_A.Lambda_curl.orig'*mtx_A.I.hat,'fro');
        ErrValidation.lemma4(7) = norm(mtx_A.Pi.tilde'*mtx_A.I.hat,'fro');
        %% Validation of Lemma 5. orthononality of Q, Q_r, P, and P_r
        x = rand(Dim.Ne,1) + 1i*rand(Dim.Ne,1);  ErrValidation.orth_Q  = norm(fnc.Qs(fnc.Q(x)) - x);
        
        x = rand(Dim.Nh,1) + 1i*rand(Dim.Nh,1);  ErrValidation.orth_P  = norm(fnc.Ps(fnc.P(x)) - x);
        x = rand(Dim.Nr,1) + 1i*rand(Dim.Nr,1);  ErrValidation.orth_Qr = norm(fnc.Qrs(fnc.Qr(x)) - x);
        x = rand(Dim.Nr,1) + 1i*rand(Dim.Nr,1);  ErrValidation.orth_Pr = norm(fnc.Prs(fnc.Pr(x)) - x);
        %% Validation of Theorem 2. SVD of C
%         x = rand(Dim.Ne,1) + 1i*rand(Dim.Ne,1); ErrValidation.SVD_C   = norm(C *x - FAME_fnc.P(Lambdas.Sigma  *fun.Qs(x))) > 1e-8;
%         x = rand(Dim.Nh,1) + 1i*rand(Dim.Nh,1); ErrValidation.SVD_Cs  = norm(C'*x - FAME_fnc.Q(Lambdas.Sigma' *fun.Ps(x))) > 1e-8;
        x = rand(Dim.Ne,1) + 1i*rand(Dim.Ne,1); ErrValidation.cSVD_C  = norm(C *x - fnc.Pr(mtx_A.Lambda_curl.Sigma.r .*fnc.Qrs(x)));
        x = rand(Dim.Nh,1) + 1i*rand(Dim.Nh,1); ErrValidation.cSVD_Cs = norm(C'*x - fnc.Qr(conj(mtx_A.Lambda_curl.Sigma.r).*fnc.Prs(x)));

    end
    
end



function y = funhand_test(b,sqrt_Sigma_r,invtildeM,tildeN)
    y = zeros(length(b),1);
    y(1:end/2) = sqrt_Sigma_r.\invtildeM(sqrt_Sigma_r.\b(1:end/2));
    
    init_y = sqrt_Sigma_r.*(tildeN(sqrt_Sigma_r.*y(1:end/2)));
    invtildeN  = @(x) pcg( tildeN, x, 1e-10, 100000, [], [], init_y/norm(init_y));

    y(end/2+1:end) = sqrt_Sigma_r.\invtildeN(sqrt_Sigma_r.\b(end/2+1:end));
end