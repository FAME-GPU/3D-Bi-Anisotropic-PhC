function fncB = fncB_Generate_Lanzcos(fnc,mtx_B,tol,iter,y)
    % 生成fncB
    % tol_M = tol(1);tol_Phi = tol(2);
    % fnc1 = @(x) [fnc.Pr(x(1:end/2) ); fnc.Qr(x(end/2+1:end))];
    % fnc2 = @(x) [mtx_B.J_int * minres(@(vec)M_int_prod_vec(mtx_B,vec),...
    %     x(1:end/2),tol_M,iter) + x(end/2+1:end); -x(1:end/2)];
    % fnc3 = @(x) [minres(@(vec)Phi_int_prod_vec(mtx_B,vec),...
    %     x(1:end/2),tol_Phi,iter); ...
    %     minres(@(vec)M_int_prod_vec(mtx_B,vec),...
    %     x(end/2+1:end),tol_M,iter)];
    % fnc4 = @(x) [minres(@(vec)M_int_prod_vec(mtx_B,vec),...
    %     mtx_B.L_int * x(1:end/2),tol_M,iter)- x(end/2+1:end); x(1:end/2)];
    % fnc5 = @(x) [fnc.Prs(x(1:end/2)); fnc.Qrs(x(end/2+1:end))];
    % fncB = fnc5( fnc4( fnc3( fnc2( fnc1( y ) ) ) ) );
    mtx_B.Phi_int =  mtx_B.N_int;% - mtx_B.J_int * (mtx_B.M_int \ mtx_B.L_int);
    fnc1 = @(x) [fnc.Pr(x(1:end/2) ); fnc.Qr(x(end/2+1:end))];
    fnc2 = @(x) [x(end/2+1:end); -x(1:end/2)];
    fnc3 = @(x) [mtx_B.Phi_int \ x(1:end/2); mtx_B.M_int \ x(end/2+1:end)];
    fnc4 = @(x) [-x(end/2+1:end); x(1:end/2)];
    fnc5 = @(x) [fnc.Prs(x(1:end/2)); fnc.Qrs(x(end/2+1:end))];    
    fncB = fnc5( fnc4( fnc3( fnc2( fnc1( y ) ) ) ) );
end

% function y = Phi_int_prod_vec(mtx_B,vec)
%     y = mtx_B.N_int * vec - mtx_B.J_int * (mtx_B.M_int \ (mtx_B.L_int * vec));
% end
function y = Phi_int_prod_vec(mtx_B,vec)
    y = mtx_B.N_int * vec - mtx_B.J_int * minres1(@(vector)M_int_prod_vec(mtx_B,vector),...
        mtx_B.L_int * vec,1e-12,1000);
end
function y = M_int_prod_vec(mtx_B,vec)
    y = mtx_B.M_int * vec;
end

    % mtx_B.Phi_int =  mtx_B.N_int - mtx_B.J_int * (mtx_B.M_int \ mtx_B.L_int);
    % fnc1 = @(x) [fnc.Pr(x(1:end/2) ); fnc.Qr(x(end/2+1:end))];
    % fnc2 = @(x) [mtx_B.J_int * (mtx_B.M_int \ x(1:end/2)) + x(end/2+1:end); -x(1:end/2)];
    % fnc3 = @(x) [mtx_B.Phi_int \ x(1:end/2); mtx_B.M_int \ x(end/2+1:end)];
    % fnc4 = @(x) [mtx_B.M_int \ (mtx_B.L_int*x(1:end/2))- x(end/2+1:end); x(1:end/2)];
    % fnc5 = @(x) [fnc.Prs(x(1:end/2)); fnc.Qrs(x(end/2+1:end))];    
    % fncB = fnc5( fnc4( fnc3( fnc2( fnc1( y ) ) ) ) );