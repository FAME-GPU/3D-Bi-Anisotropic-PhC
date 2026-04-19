function [ vec_t ] = SIRA_solve_residual_LS_GEP(A, B, sigma, rhs, LS_solver, pcg_iter_number, ...
    linear_system_tol,no_file, prec_M, init_vec)
global time_minres flag_gpu
%Solve the linear system
%        ( A - tau I ) vec_t = rhs.
%    
%49行开始输出   改不对称时，添加55-59行gk
global total_inner_iter 

if ( nargin <= 7 )
    switch LS_solver
        case 'pcg'
            [vec_t, flag, RELRES, iter_num ] = pcg(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
    
        case 'minres'
            [vec_t, flag, RELRES, iter_num ] = minres(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'bicgstab'
            [vec_t, flag, RELRES, iter_num ] = bicgstab(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'bicgstabl'
            [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'cgs'
            [vec_t, flag, RELRES, iter_num ] = cgs(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'gmres'
            RESTRT       = 10;
            [vec_t, flag, RELRES, iter_num2 ] = gmres(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                RESTRT, linear_system_tol, pcg_iter_number); 
            iter_num = iter_num2(1)*RESTRT+iter_num2(2);
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num2(1)*RESTRT+iter_num2(2), RELRES, flag);
            
        otherwise
            [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
    end
elseif (nargin == 8)
    switch LS_solver
        case 'minres'
            minres_time = tic;            
            [vec_t, flag, RELRES, iter_num ] = minres1(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number); 
            if flag_gpu
                wait(gpuDevice);
            end
            time_minres.minres = time_minres.minres + toc(minres_time);
            fprintf(no_file,'iter_minres = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        case 'gmres'
             RESTRT = 5;
            [vec_t, flag, RELRES, iter_num2 ] = gmres(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                RESTRT, linear_system_tol, pcg_iter_number); 
            iter_num = iter_num2(1)*RESTRT+iter_num2(2);
            fprintf(no_file,'iter_gmres = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
    end
elseif ( nargin == 9 )
    switch LS_solver
        case 'pcg'
            [vec_t, flag, RELRES, iter_num ] = pcg(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, prec_M); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
    
        case 'minres'
            %maxit = min(250,pcg_iter_number);
            [vec_t, flag, RELRES, iter_num ] = minres(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, prec_M); 
            fprintf(no_file,'iter_minres = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'bicgstab'
            [vec_t, flag, RELRES, iter_num ] = bicgstab(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, prec_M); 
            fprintf('iter_bicgstab = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'bicgstabl'
            [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, prec_M); 
            fprintf('iter_bicgstabl = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'cgs'
            [vec_t, flag, RELRES, iter_num ] = cgs(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, prec_M); 
            fprintf('iter_cgs = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
            
            if ( flag ~= 0 )
                RESTRT       = 30; 
                [vec_t, flag, RELRES, iter_num ] = gmres(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                    RESTRT, linear_system_tol, 5, prec_M, [], vec_t); 
                fprintf('iter_gmres = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num(1)*RESTRT+iter_num(2), RELRES, flag);
            end
            
        case 'tfqmr'
            [vec_t, flag, RELRES, iter_num ] = tfqmr(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, prec_M); 
            fprintf('iter_tfqmr = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'gmres'
            RESTRT       = 30;
            maxiter      = min(6,round(pcg_iter_number / RESTRT));
            [vec_t, flag, RELRES, iter_num2 ] = gmres(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                RESTRT, linear_system_tol, maxiter, prec_M); 
            iter_num = iter_num2(1)*RESTRT+iter_num2(2);
            fprintf('iter_gmres = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num2(1)*RESTRT+iter_num2(2), RELRES, flag);
        
        otherwise
            [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, prec_M); 
            fprintf('iter_bicgstabl = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
            
    end  
    
    if ( iter_num < 4 )
        [vec_t, flag, RELRES, iter_num ] = minres(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number); 
        fprintf('iter_minres = %4.0f without preconditioner; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
    end
    
else
    switch LS_solver
        case 'pcg'
            [vec_t, flag, RELRES, iter_num ] = pcg(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, [], [], init_vec); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
    
        case 'minres'
            [vec_t, flag, RELRES, iter_num ] = minres(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, [], [], init_vec); 
            fprintf(no_file,'iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'bicgstab'
            [vec_t, flag, RELRES, iter_num ] = bicgstab(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, [], [], init_vec); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'bicgstabl'
            [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, [], [], init_vec); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'cgs'
            [vec_t, flag, RELRES, iter_num ] = cgs(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, [], [], init_vec); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'gmres'
            RESTRT       = 10;
            [vec_t, flag, RELRES, iter_num2 ] = gmres(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                RESTRT, linear_system_tol, pcg_iter_number, [], [], init_vec); 
            iter_num = iter_num2(1)*RESTRT+iter_num2(2);
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num2(1)*RESTRT+iter_num2(2), RELRES, flag);
            
        otherwise
            [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, [], [], init_vec); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
    end
    
end

total_inner_iter = total_inner_iter + iter_num;
%outer_iter_count = outer_iter_count + 1;

end

% =========================================================================
%
% function sol = Precond_M_LS_GEP(x, sigma, prec_M)
% %UNTITLED2 Summary of this function goes here
% %   Detailed explanation goes here
% 
% sol = prec_M(x, sigma);
% 
% end
