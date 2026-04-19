function y = GEP_shift_mtx_prod_vec(vec, A, B, sigma )
global time_minres flag_gpu
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

mtx_x_time = tic;
% lambda = 1e-3;
% y = MtxVecMult_A(vec) - sigma * MtxVecMult_B(vec) + lambda * vec;
y = MtxVecMult_A(vec);
%y = MtxVecMult_A(vec) - sigma * MtxVecMult_B(vec);
if flag_gpu
    wait(gpuDevice);
end
time_minres.mtx_prod_vec.all = time_minres.mtx_prod_vec.all + toc(mtx_x_time);

end

% =======================================
%
 function u = Amtrix(x, A)
   u = A * x;
 end
 
 function u = Bmtrix(x, B)
   u = B * x;
 end
