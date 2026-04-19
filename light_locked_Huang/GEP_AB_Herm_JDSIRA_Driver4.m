function [EW, EV] = GEP_AB_Herm_JDSIRA_Driver4 (A, B, mtxdim, no_restart, ew_number, ...
                stop_tolerance, initial_vec, target, no_file, target_type, CorEq, ...
                LSinfo, Precond_M_LS, Construct_precond_LU_JD, opt, mtx_B, fnc, wave_vec, N, gamma, model_type)

% mtx_A      : symmetric cofficient matrix of the eigenproblem                               
% mtxdim     : dimension of the cofficient matricx mtx_A      
% no_restart : number to restart the JD iteration        
% ew_number  : total number of eigenpairs wanted        
global zflag_ew flag_no_ew dim_eigsp total_inner_iter
global flag_gpu time_minres

    RPPD                      = max(4, size(initial_vec,2)); 
    ProjProbDim               = size(initial_vec,2);
    % init_Vec                  = zeros(mtxdim,ew_number+RPPD-1);
    init_Vec                  = zeros(mtxdim,ew_number+RPPD-1,'double', 'gpuArray');
    % if flag_gpu
    %     init_Vec = gpuArray(init_Vec);
    % end
    init_Vec(:,1:ProjProbDim) = initial_vec; 
    
    flag_no_ew    = 0;
    dim_eigsp     = 0;
    % EV = zeros(mtxdim,ew_number); 
    EV = zeros(mtxdim,ew_number,'double', 'gpuArray');
    % if flag_gpu
    %     EV = gpuArray(EV);
    % end

    RestartProjProbDim = RPPD;       % RestartProjProbDim must >= 1  

    % zflag_ew           = zeros(ew_number+RPPD,1);
    zflag_ew           = zeros(ew_number+RPPD,1,'double', 'gpuArray');
    % if flag_gpu
    %     zflag_ew = gpuArray(zflag_ew);
    % end

    % --- iteration for finding the first ew_number eigenpairs

    no_computed_ew = 1; 

    stop_flag = 1;
    
    while ( no_computed_ew <= ew_number && stop_flag ~= 0 ) 

          % --- compute the no_ite_ew-th eigenpair

          total_inner_iter = 0;
          compute_an_eig_time = tic;  
          [no_computed_ew, EV, init_Vec, ProjProbDim, stop_flag] = GEP_AB_Herm_JDSIRA_Locking4 ...
              ( A, B, EV, no_computed_ew, init_Vec, ProjProbDim, target, target_type, ...
              no_file, RestartProjProbDim, no_restart, stop_tolerance, ew_number, CorEq, ...
              LSinfo, Precond_M_LS, Construct_precond_LU_JD);
          if flag_gpu
             wait(gpuDevice);
          end
          time_minres.compute_an_eig = time_minres.compute_an_eig + toc(compute_an_eig_time);

           if ( ProjProbDim == 0 ) 
               init_Vec(:,1) = randn(size(init_Vec,1),1); 
               rsdl          = norm(init_Vec(:,1));
               init_Vec(:,1) = init_Vec(:,1) / dsqrt(rsdl); 
               ProjProbDim   = 1;
           end
           

           fprintf(no_file,'total_inner_iter(kk,%2.0f) = %7.0f; \n', no_computed_ew-1, total_inner_iter);
           fprintf('%d eigenvalues have been caculated.\n', no_computed_ew-1);
   
           result.opt = opt;
           ew = zflag_ew(1:no_computed_ew-1,1);
           [ result.Eigenmode.e1, result.Eigenmode.e2, result.Eigenmode.e3, ...
            result.Eigenmode.h1, result.Eigenmode.h2, result.Eigenmode.h3] ...
             = Eigen_Restoration_Chiral_no_inv(EV(:,1:no_computed_ew-1), opt.comp.Dim, mtx_B, fnc.Qr, fnc.Pr);
           [result.res_NFSEP, result.res_ORGSEP] = EvaluteError_Chiral_no_inv(EV(:,1:no_computed_ew-1), ew,result.Eigenmode, opt, B, A, mtx_B);
           filename = ['./test_data_store/' model_type '_n' num2str(N) '_gamma_choose4_=' sprintf('%.7f', gamma) '.mat'];
            save(filename, 'ew', 'result', 'wave_vec', 'time_minres');

    end 
    
    EW = zflag_ew(1:ew_number,1);
    fprintf(no_file,'kk = kk +1; \n\n');
end