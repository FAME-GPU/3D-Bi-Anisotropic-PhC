function [ e1, e2, e3, h1, h2, h3] = Eigen_Restoration_Chiral_Lanzcos( EV_NF, Dim, mtx_B, fnc_Qr, fnc_Pr)
%  利用NFSEP的特征向量EV_NF生成ORGSEP的特征向量
%   function [ e1, e2, e3, h1, h2, h3] = Eigen_Restoration_Chiral( EV, lambda, Dim, mtx_B, fnc_Qr, fnc_Pr, Sigma_r)
    tol = [1e-15;1e-13];iter = 5000;
    fncB_1 = @(x)fncB1_Generate(mtx_B,tol,iter,x);
    fncB_2 = @(x)fncB2_Generate(mtx_B,tol,iter,x);
    fncB_inv = @(x) -1i * fncB_2( fncB_1(x) );
	fnc_tmp1 = @(x) [fnc_Pr( x(1:end/2) ); fnc_Qr( x(end/2+1:end) )];
	fnc_invariant = @(x) fncB_inv( fnc_tmp1(x) );
    
    %% Restore eigenpairs
    omega = zeros(size(EV_NF,2),1);
    for j = 1 : size(EV_NF, 2)
        v = EV_NF(:,j);
        % b = 1i * [ -fnc_Qr( Sigma_r .* v(end/2+1:end) ); fnc_Pr( Sigma_r.*v(1:end/2) ) ] / lambda(j);
% 		e_h_j = fnc_invariant( sigma_r_sqrt_double .\ v );
        e_h_j = fnc_invariant( v );
        
        % [e_h_j, flag, relres, iter]= minres(fnc_handle, b, 1e-6, LS_maxit);

        e = e_h_j(1:end/2);
        h = e_h_j(end/2+1:end);
        
        e1(:,j) = e(                1:Dim.Ne1        );
        e2(:,j) = e(        Dim.Ne1+1:Dim.Ne1+Dim.Ne2);
        e3(:,j) = e(Dim.Ne1+Dim.Ne2+1:end            );
        h1(:,j) = h(                1:Dim.Nh1        );
        h2(:,j) = h(        Dim.Nh1+1:Dim.Nh1+Dim.Nh2);
        h3(:,j) = h(Dim.Nh1+Dim.Nh2+1:end            );
    end
end