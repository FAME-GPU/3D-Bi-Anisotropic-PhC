function [ew_order, EV, init_Vec, ProjProbDim, stop_flag] = GEP_AB_Herm_JDSIRA_Locking4(A, B, EV, ew_order, ...
    init_Vec, ProjProbDim, target, Choose_Type, no_file, RestartProjProbDim, no_restart, ...
    vec_r_tol, ew_number, CorEq_org, LSinfo, Precond_M_LS, Construct_precond_LU_JD )
%
% Use Jacobi-Davidson (JD) or Shift-Invert Residual Arnoldi (SIRA) Method to 
% compute the target eigenpair of the generalized eigenvalue problem
%                A x = lambda B x
% where A and B are Hermitian.
% 181、306行输出
global theta ev_crnt 
global zflag_ew flag_no_ew dim_eigsp
global time_minres
global flag_gpu

if isa(A,'function_handle')
    MtxVecMult_A = A;
else
    % Check matrix and right hand side vector inputs have appropriate sizes
    [m,n] = size(A);
    if (m ~= n)
        fprintf('error:NonSquareMatrix\n');return;
    end
    MtxVecMult_A = @(x)Amtrix(x,A); 
end

if isa(B,'function_handle')
    MtxVecMult_B = B;
else
    % Check matrix and right hand side vector inputs have appropriate sizes
    [m,n] = size(B);
    if (m ~= n)
        fprintf('error:NonSquareMatrix\n');return;
    end
    MtxVecMult_B = @(x)Bmtrix(x,B); 
end

    % -------------------------
    % initialize matrix V, W, M
    % -------------------------
    
    if ( ~isfield(LSinfo,'maxit') ) 
        max_iter_no = 500;
    else
        max_iter_no = LSinfo.maxit;
    end
    
    if ( ~isfield(LSinfo,'tol') )
        linear_system_tol = 1.0e-2;
    else
        linear_system_tol = LSinfo.tol;
    end
    
    tic;
    stop_flag          = 0; 
    EPIteNo            = 1; 
    iter_flag          = 0;
    iter               = 60; %15;
    pcg_iter_number    = max_iter_no; %1000;
    prev_rsdl          = 1.0e2; 
    no_ew_near_tol     = 0;
    mtxdim             = size(init_Vec,1);
    target_ritz_ew     = target;
    change_target      = 0;
    prev_ritz_ew       = 1.0e14;
    tol_ritz_ew        = 5.0e-3;
    Choose_Type_org    = Choose_Type;

    mtxV                  = zeros( mtxdim,       no_restart+2 ,'double','gpuArray');
    mtxM                  = zeros( no_restart+2, no_restart+2 ,'double','gpuArray');
    mtxL                  = zeros( no_restart+2, no_restart+2 ,'double','gpuArray');
    % mtxV                  = zeros( mtxdim,       no_restart+2 ); 
    % mtxM                  = zeros( no_restart+2, no_restart+2 );
    % mtxL                  = zeros( no_restart+2, no_restart+2 );
    % if flag_gpu
    %     mtxV = gpuArray(mtxV);
    %     mtxM = gpuArray(mtxM);
    %     mtxL = gpuArray(mtxL);
    % end
%     init_Vec(:,1) = init_Vec(:,1) / sqrt(init_Vec(:,1)' * MtxVecMult_B(init_Vec(:,1)));
    mtxV(:,1:ProjProbDim) = init_Vec(:,1:ProjProbDim);
    
    while ( EPIteNo <= 6000 && stop_flag == 0 )

       %
       % Compute the Rayleigh quotient
       %      V^{*} A V   and   V^{*} B V
       %
       % MtxVecMult( V(:,ii), 0 ) : Compute A * V(:,ii)
       % MtxVecMult( V(:,ii), 1 ) : Compute B * V(:,ii)
       %
       for ii = 1:ProjProbDim
           mtx_A_prod_V              = MtxVecMult_A( mtxV(:,ii) );  %????????????????
           mtxM(ii:ProjProbDim,ii)   = mtxV(:,ii:ProjProbDim)' * mtx_A_prod_V;
           mtxM(ii,ii+1:ProjProbDim) = mtxM(ii+1:ProjProbDim,ii)';
           
           mtx_A_prod_V              = MtxVecMult_B( mtxV(:,ii) );  
           mtxL(ii:ProjProbDim,ii)   = mtxV(:,ii:ProjProbDim)' * mtx_A_prod_V;
           mtxL(ii,ii+1:ProjProbDim) = mtxL(ii+1:ProjProbDim,ii)';
           % if mtxL(ii:ProjProbDim,ii) == 0                           %gk correct
           %     mtxL(ii:ProjProbDim,ii) = mtxL(ii:ProjProbDim,ii);%+1e-12*1i;
           %     mtxL(ii,ii+1:ProjProbDim) = mtxL(ii+1:ProjProbDim,ii)';
           % end
       end 

       % -----------------------------
       % iterate to find the eigenpair
       % -----------------------------

       iteno      = ProjProbDim; 
       while (iteno <= no_restart && stop_flag == 0 ) 

          % --- solve for the expanded generalized eigensystem
          time_eig = tic;
          % [VRR, ew] = eig(mtxM(1:iteno, 1:iteno), mtxL(1:iteno, 1:iteno));
          [VRR, ew] = eig(gather(mtxM(1:iteno, 1:iteno)), gather(mtxL(1:iteno, 1:iteno)));
          if flag_gpu
              wait(gpuDevice);
          end
          time_minres.sira.eig = time_minres.sira.eig + toc(time_eig);
          ew        = diag(ew);
          
          if flag_gpu
              ew = gpuArray(ew);
              VRR = gpuArray(VRR);
          end
          
          if (sum(isinf(ew))>1)
              error("inf merg.")
          end
          
%           if sum(isinf(ew)) > 0
%               break;
%           end
          % --- select the desired eigenpair (lambda,x) 
          %     where lambda is the
          %     closest to (real_target, imag_target)
%           if strcmp(Choose_Type, 'RGTC')
%             test_ew = ew;
%             test_ew(real(test_ew) < real(target_ritz_ew)) = 1e3;
%             [~,idx] = sort(real(test_ew-target_ritz_ew));
%             idx_order = idx(ew_order);
%             ritz_ew = ew(idx);
%             ritz_ev = VRR(:,idx);
%             ewno = ew_order;
%           elseif strcmp(Choose_Type, 'CL')
%               test_ew = ew;
% %             test_ew(real(test_ew) < real(target_ritz_ew)) = 1e3;
%             [~,idx] = sort(abs(test_ew-target_ritz_ew));
%             idx_order = idx(ew_order);
%             ritz_ew = ew(idx);
%             ritz_ev = VRR(:,idx);
%             ewno = ew_order;
%           end
          
          time_Choose_Ritz_Pair_Unsymm = tic;
          [ewno, ritz_ew, ritz_ev] = Choose_Ritz_Pair_Unsymm4(iteno, ew, ew_order, ...
              target_ritz_ew, VRR, vec_r_tol, Choose_Type, RestartProjProbDim, ...
              target, Choose_Type_org);  
          if flag_gpu
            wait(gpuDevice);
          end
          time_minres.sira.Choose_Ritz_Pair = time_minres.sira.Choose_Ritz_Pair + toc(time_Choose_Ritz_Pair_Unsymm);
                  
          % ----------------------------------------------------
          % comput u (eigenvector of the PEP), p = A'(labda)u, 
          %        r = A(lambda)u (residual of the PEP)
          % ----------------------------------------------------

          %         
          %  ev_crnt : current target approximated eigenvector
          %
          ev_crnt  = mtxV(:,1:iteno) * ritz_ev(1:iteno,1);  
          ev_crnt  = ev_crnt / norm( ev_crnt); 
          %
	  % Compute the residual of the ritz pair ( ritz_ew(1), ev_crnt )
	  %
          rsdl_vec = MtxVecMult_A(ev_crnt);
          rsdl_vec = rsdl_vec - ritz_ew(1) * MtxVecMult_B(ev_crnt); 
          rsdl     = norm(rsdl_vec);
  
          print_result = 1;
          if ( print_result == 1 )
             if ( imag(ritz_ew(1)) > 0.0  ) 
                fprintf(no_file,'ew_crnt(%2.0f,%4.0f) = %18.12e+%18.12ei; rsdl(%2.0f,%4.0f) = %8.2e;\n', ...
                   ew_order, EPIteNo, real(ritz_ew(1)), imag(ritz_ew(1)), ...
                   ew_order, EPIteNo, rsdl);
             elseif ( imag(ritz_ew(1)) < 0  )
                fprintf(no_file,'ew_crnt(%2.0f,%4.0f) = %18.12e-%18.12ei; rsdl(%2.0f,%4.0f) = %8.2e;\n', ...
                   ew_order,EPIteNo,real(ritz_ew(1)),abs(imag(ritz_ew(1))), ...
                   ew_order, EPIteNo, rsdl);
             else
                fprintf(no_file,'ew_crnt(%2.0f,%4.0f) = %18.12e; rsdl(%2.0f,%4.0f) = %8.2e;\n', ...
                   ew_order, EPIteNo, real(ritz_ew(1)), ...
                   ew_order, EPIteNo, rsdl);
             end 
          end

          if ( rsdl <= 10.0e0 * vec_r_tol && rsdl > vec_r_tol )  
             no_ew_near_tol = no_ew_near_tol + 1;
          else
             no_ew_near_tol = 0;
          end

          % -----------------
          % check convergence
          % -----------------

          if ( rsdl <= vec_r_tol || no_ew_near_tol >= 5 )    % 残差连续多次接近阈值，标记收敛
             stop_flag = 1;     
          else  
             if ( rsdl < tol_ritz_ew && change_target == 1 && ...
                  abs(prev_ritz_ew - ritz_ew(1))/abs(ritz_ew(1)) < 5.0e-2 )

                target_ritz_ew = ritz_ew(1);
                change_target  = 0;
                Choose_Type    = 'CL';
             end
             prev_ritz_ew = ritz_ew(1);
%
%            modify at 2007, 1, 19
%
%             if ( iteno + increment <= no_restart ) then
             if ( iteno + 1 <= dim_eigsp+no_restart )  
%
%               Compute the appending vector
%
                %
                % Adaptive determine the iteration numbers for solving correction 
                % equation in JD method
                %

%                 if ( rsdl >= 1.0e-1 && EPIteNo > 14 ) 
%                    iter            = 30;
%                    iter_flag       = 1;
%                    if ( pcg_iter_number < 1500 )
%                        pcg_iter_number = 1500;
%                    end
%                 elseif ( iter_flag == 1 )  
%                    iter            = 60; %15;
%                    iter_flag       = 0;
%                    if ( pcg_iter_number < 1000 )
%                        pcg_iter_number = 1000;
%                    end
%                 elseif ( rsdl < 1.0e-1 && prev_rsdl/rsdl < 4.0e0 ) 
%                    iter            = min(300, iter+10); %min(200, iter+5);
%                    pcg_iter_number = min(2000,pcg_iter_number+100);Error for multiplicity of convergent
%                 end 
                
                if ( rsdl < 1.0e-1 && prev_rsdl/rsdl < 4.0e0 )    % rsdl<0.1时调整迭代次数
                   iter            = min(300, iter+10); %min(200, iter+5);
                   pcg_iter_number = min(floor(max_iter_no*1.2),pcg_iter_number+10);
                end 
                
                %pcg_iter_number = 1500;
                prev_rsdl = rsdl;

%
%                  init_Vec(:,1) := M_A^{-1} r = M_A^{-1} tmpVt(:,2)
%                  tmpVt(:,2)    := M_A^{-1} p = M_A^{-1} tmpVt(:,1)
%
                theta     = ritz_ew(1);
                LS_solver = LSinfo.solver;
                
                if ( ew_order == 1 && iteno <= 10 && rsdl > 1.0e-1 )
                    CorEq = 'SIRA';
                    if ( strcmp(CorEq_org,'JD') || strcmp(CorEq_org,'jd') )
                        %LSinfo.precond = 'yes';
                        precond_info    = 'no'; %'yes';
                        pcg_iter_number = 500;
                    else
                        precond_info    = 'no';
                        pcg_iter_number = 5000;
                        LS_solver       = 'minres';
                    end
                    iter_flag = 0;
%                     LSinfo.precond = 'yes';
%                     CorEq          = CorEq_org;
%                     iter           = 30;
                else
                    if ( iter_flag == 0 )
                        pcg_iter_number = max_iter_no;
                        iter_flag       = 1;
                    end
                    CorEq        = CorEq_org; %'JD'; %'SIRA';
                    precond_info = LSinfo.precond;
                end
                %CorEq = CorEq_org;
                
                switch CorEq
                    case 'JD' 
                        LU_precond    = Construct_precond_LU_JD(theta);
                        LS_solver     = LSinfo.solver;
                        init_Vec(:,1) = solve_approx_vec_t_GEP_Herm(A, B, @(x)solve_Minv_b( x, LU_precond ), LS_solver, theta, ...
                            rsdl_vec, ev_crnt, ew_order, pcg_iter_number, rsdl);
%                         init_Vec(:,1) = solve_approx_vec_t_GEP_Herm_test(A, B, Construct_precond_LU_JD, LS_solver, theta, ...
%                             rsdl_vec, ev_crnt, ew_order, iter, rsdl); 
                           
                    case 'SIRA' 

%                         linear_system_tol = 1.0e-3;
                        %pcg_iter_number   = 2*10^3; %iter; %2 * 10^3;
%                         LS_solver         = LSinfo.solver; %'minres';

%                         save rhs_LS_dat rsdl_vec target

                        if ( strcmp(precond_info, 'yes') || strcmp(precond_info, 'YES') ) 
                            init_Vec(:,1) = SIRA_solve_residual_LS_GEP(A, B, target, rsdl_vec, ...
                                LS_solver, pcg_iter_number, linear_system_tol, Precond_M_LS,no_file);
                        %        LS_solver, pcg_iter_number, linear_system_tol, ev_crnt, Precond_M_LS); 
                        else 
                            init_Vec(:,1) = SIRA_solve_residual_LS_GEP(A, B, target, rsdl_vec, ...
                                LS_solver, pcg_iter_number, linear_system_tol,no_file); 
%                                LS_solver, pcg_iter_number, linear_system_tol, ev_crnt);
                        end
                end

                rsdl          = norm(init_Vec(:,1));
                if ( rsdl < 1.0e-9 )  
                   init_Vec(:,1) = init_Vec(:,1) / rsdl;
                end


                % ----------------------------------------------
                % orthogonalize t against V and then nomalize it
                % ----------------------------------------------
                for j = 1:iteno
                   ew_tmp        = mtxV(:,j)' * init_Vec(:,1);
                   init_Vec(:,1) = init_Vec(:,1) - ew_tmp * mtxV(:,j);
                end  
                mtxV(:,iteno+1) = init_Vec(:,1) / norm(init_Vec(:,1));
 
                % ---------------
                % update V^T A V  and  V^T B V
                % ---------------
                mtx_A_prod_V            = MtxVecMult_A( mtxV(:,iteno+1) ); 
                mtxM(1:iteno+1,iteno+1) = mtxV(:,1:iteno+1)' * mtx_A_prod_V;
                mtxM(iteno+1,1:iteno)   = mtxM(1:iteno,iteno+1)';
                
                mtx_A_prod_V            = MtxVecMult_B( mtxV(:,iteno+1) ); 
                mtxL(1:iteno+1,iteno+1) = mtxV(:,1:iteno+1)' * mtx_A_prod_V;
                mtxL(iteno+1,1:iteno)   = mtxL(1:iteno,iteno+1)';
             end 

          end 

          iteno   = iteno + 1;
          EPIteNo = EPIteNo + 1;
          
       end

       % --- construct the initial search space with ProjProbDim Ritz pairs 
       %     for restarting
       %
       % The orthonormal basis for the convergent eigenvactors is always locked at 
       % the first "dim_eigsp" columns of the matrix mtxV. The current approximate 
       % eigenvector is put on the "dim_eigsp"+1-th column of mtxV. 
       % We also collect the Ritz vectors for the Ritz values which are the second, 
       % third, ... closest to the target value to accelerate the convergence in computing	
       % next target eigenpair.
       %
       if ( stop_flag == 0 )

          mtxM  = zeros(no_restart+2, no_restart+2,'double','gpuArray'); 
          % if flag_gpu
          %     mtxM = gpuArray(mtxM);
          % end
          iteno = iteno - 1;

          jj               = min(RestartProjProbDim-1,ew_number+RestartProjProbDim-1);
          init_Vec(:,1:jj) = mtxV(:,1:iteno) * ritz_ev(1:iteno,2:jj+1);

          k             = dim_eigsp; 
          mtxV(:, k+1)  = ev_crnt;
          k             = k + 1; 

          mtxV(:,k+1:k+jj) = init_Vec(:,1:jj);
          k                = k + jj;

%           i            = dim_eigsp; 
          [i, mtxV(:,1:i)] = Gram_Schmidt(k, mtxV(:, 1:k));
          ProjProbDim      = i;

       end 
    end 

    %
    %   Compute Ritz pairs and the corresponding residual
    %
    rsdl_check = 1;
    tmp_rsdl   = ones(no_restart+2, no_restart+2,'double','gpuArray');
    % if flag_gpu
    %     tmp_rsdl = gpuArray(tmp_rsdl);
    % end

    jj         = 0;
    for i = 2:ewno %RestartProjProbDim
          if ( jj+1 <= ew_number+RestartProjProbDim-1 )   
             init_Vec(:, jj+1) = mtxV(:,1:iteno-1) * ritz_ev(1:iteno-1,i); 
             init_Vec(:, jj+1) = init_Vec(:,jj+1) / norm(init_Vec(:,jj+1)); 
             jj                = jj + 1;
             if ( rsdl_check == 1 ) 
                 rsdl_vec   = MtxVecMult_A(init_Vec(:,jj)) - ritz_ew(i) * MtxVecMult_B(init_Vec(:,jj));
                tmp_rsdl(i) = norm(rsdl_vec); 
                
                if ( imag(ritz_ew(i)) > 0 ) 
                   fprintf(no_file,'rsdl(%24.16e+%24.16ei) = %14.4e \n', real(ritz_ew(i)), ...
                       imag(ritz_ew(i)),tmp_rsdl(i));
                elseif ( imag(ritz_ew(i)) < 0 )
                   fprintf(no_file,'rsdl(%24.16e-%24.16ei) = %14.4e \n', real(ritz_ew(i)), ...
                       -imag(ritz_ew(i)),tmp_rsdl(i));
                else
                   fprintf(no_file,'rsdl(%24.16e) = %14.4e \n', real(ritz_ew(i)),tmp_rsdl(i));
                end
             end
          end
    end
    
    flag_no_ew = flag_no_ew + 1;
    k          = dim_eigsp;

    %
    % Save the converging eigenpair and lock it on matrix V
    %

       tmpVt(:,1)           = ev_crnt / norm(ev_crnt);
       tmpVt2 = zeros(size(ev_crnt,1),1,'double','gpuArray');
       % if flag_gpu
       %     tmpVt2 = gpuArray(tmpVt2);
       % end
       tmpVt(:,2)           = tmpVt2;
       if ( rsdl_check == 1 ) 
           rsdl_vec = MtxVecMult_A(tmpVt(:,1)) - ritz_ew(1) * MtxVecMult_B(tmpVt(:,1));
           rsdl     = norm(rsdl_vec);
           
           if ( imag(ritz_ew(1)) > 0 ) 
               fprintf(no_file,'rsdl(%24.16e+%24.16ei) = %14.4e \n', real(ritz_ew(1)), ...
                       imag(ritz_ew(1)),rsdl);
           elseif ( imag(ritz_ew(1)) < 0 )
               fprintf(no_file,'rsdl(%24.16e-%24.16ei) = %14.4e \n', real(ritz_ew(1)), ...
                       -imag(ritz_ew(1)),rsdl);
           else
               fprintf(no_file,'rsdl(%24.16e) = %14.4e \n', real(ritz_ew(1)),rsdl);
           end
                
       end 
       mtxV(:,k+1  )        = tmpVt(:,1);
       EV(:,k+1)            = tmpVt(:,1);
       k                    = k + 1;
       zflag_ew(flag_no_ew) = ritz_ew(1);
       dim_eigsp            = dim_eigsp + 1;
    %
    %   Check any other converging eigenpair. If yes, lock it again into V
    %
    % flag_cn = zeros(ewno,1);
    flag_cn = zeros(ewno,1,'double','gpuArray');
    % if flag_gpu
    %     flag_cn = gpuArray(flag_cn);
    % end
    j1      = 0;
    j       = 0;
    for i = 2:ewno %RestartProjProbDim 
          if ( j1+1 <= ew_number+RestartProjProbDim-1 ) 
             if ( tmp_rsdl(i) <= vec_r_tol)
                mtxV(:,k+1)          = init_Vec(:,j1+1);
                EV(:,k+1)            = init_Vec(:,j1+1);
                k                    = k + 1;
                flag_cn(i)           = 1;
                j                    = j + 1;
                flag_no_ew           = flag_no_ew + 1;
                zflag_ew(flag_no_ew) = ritz_ew(i);
                dim_eigsp            = dim_eigsp + 1;
                j1                   = j1 + 1;
                output(no_file, ew_order+j, 0, 0.0d0, ritz_ew(i), tmp_rsdl(i));
             end
          end 
    end

    % --- construct the initial search space
    %
    % Now, the first k columns of mtxV are the convergent eigenvectors.
    % The good initial vectors will be appended into mtxV from k+1-th column.
    %
    j1 = 0;
    for i = 2:ewno 
          if ( flag_cn(i) == 0 )
             if ( j1+1 <= ew_number+RestartProjProbDim-1 ) 
                mtxV(:,k+1) = init_Vec(:,j1+1);
                k           = k + 1;
             end
          end
          j1 = j1 + 1; 
    end

    fprintf('\n');
    
    %
    % Find the orthonormal basis of mtxV and reset it into init_vec
    %
    [i, mtxV(:,1:i)]          = Gram_Schmidt(k, mtxV(:, 1:k));
    ProjProbDim               = min(i, ew_number+RestartProjProbDim-1);
    init_Vec(:,1:ProjProbDim) = mtxV(:,1:ProjProbDim);

    output(no_file, ew_order, EPIteNo, toc, ritz_ew(1), rsdl);

    ew_order                  = ew_order + j + 1; 

end

% =======================================
%
 function u = Amtrix(x, A)
   u = A * x;
 end
 
 function u = Bmtrix(x, B)
   u = B * x;
 end

