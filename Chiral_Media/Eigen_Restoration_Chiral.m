function [ e1, e2, e3, h1, h2, h3] = Eigen_Restoration_Chiral( EV, Dim, mtx_B, fnc_Qr, fnc_Pr)
%   function [ e1, e2, e3, h1, h2, h3] = Eigen_Restoration_Chiral( EV, lambda, Dim, mtx_B, fnc_Qr, fnc_Pr, Sigma_r)
    L = mtx_B.L_int; % zeta
    J = mtx_B.J_int; % xi
	M_inv = mtx_B.M_inv_int;
	Phi_inv = mtx_B.Phi_inv_int;
	
	fncB_1 = @(x) [Phi_inv * (J * (M_inv * x(1:end/2))) + Phi_inv * x(end/2+1:end); -M_inv * x(1:end/2)];
	fncB_2 = @(x) [x(1:end/2); -M_inv*(L * x(1:end/2)) + x(end/2+1:end)];
	fncB_inv = @(x) -1i * fncB_2( fncB_1(x) );
	
% 	sigma_r_sqrt = sqrt(Sigma_r);
% 	fnc_tmp1 = @(x) [fnc_Pr( sigma_r_sqrt .* x(1:end/2) ); fnc_Qr( sigma_r_sqrt .* x(end/2+1:end) )];
	fnc_tmp1 = @(x) [fnc_Pr( x(1:end/2) ); fnc_Qr( x(end/2+1:end) )];
	fnc_invariant = @(x) fncB_inv( fnc_tmp1(x) );
    
% 	sigma_r_sqrt_double = kron(ones(2,1), sigma_r_sqrt);
    %% Restore eigenpairs
    omega = zeros(size(EV,2),1);
    for j = 1 : size(EV, 2)
        v = EV(:,j);
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

