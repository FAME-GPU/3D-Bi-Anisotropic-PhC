function [Lambda,Pi,I] = SVD_Matrix_curl(Lambda_1D, Lambda_3D, Dim)
    I_m       = @(m)   speye(m);
    O_mn      = @(m,n) sparse(m,n);
    
    n1 = Dim.n1; n2 = Dim.n2; n3 = Dim.n3;
    m1 = Dim.m1; m2 = Dim.m2; m3 = Dim.m3;
    
    row.x = n1; row.y = n2; row.z = n3;
    col.x = m1; col.y = m2; col.z = m3;
    
    Ne1 = Dim.Ne1; Ne2 = Dim.Ne2; Ne3 = Dim.Ne3;
    Nh1 = Dim.Nh1; Nh2 = Dim.Nh2; Nh3 = Dim.Nh3;
    
    Oxx = O_mn(Nh1,Ne1);
    Oyy = O_mn(Nh2,Ne2);
    Ozz = O_mn(Nh3,Ne3);
    
    I.Psi   = DistOp.I_Psi(row,col);
    I.Phi   = DistOp.I_Phi(row,col);
    I.hat   = DistOp.I_hat(row,col);
    I.tilde = DistOp.I_tilde(row,col);
 
    Lambda.orig = ...
        [            Oxx, -Lambda_3D.C13,  Lambda_3D.C12  ;
           Lambda_3D.C23,            Oyy, -Lambda_3D.C21  ;
          -Lambda_3D.C32,  Lambda_3D.C31,            Ozz ];
      
    O = sparse(Dim.Nq,Dim.Nq);
    Lambda.r = ...
        [            O, -Lambda_3D.z,  Lambda_3D.y ;
           Lambda_3D.z,            O, -Lambda_3D.x ;
          -Lambda_3D.y,  Lambda_3D.x,            O];  
    
    Lambda.tilde_q = blkdiag(...
         kron( Lambda_1D.z(Dim.rmc_1D.z+1:end,1:end)'*Lambda_1D.z(Dim.rmc_1D.z+1:end,1:end) , kron( I_m(m2) , I_m(Dim.rmc_1D.x) ) ) +...
         kron( I_m(m3) , kron( Lambda_1D.y(Dim.rmc_1D.y+1:end,1:end)'*Lambda_1D.y(Dim.rmc_1D.y+1:end,1:end) , I_m(Dim.rmc_1D.x) ) ) ,... 
         kron( Lambda_1D.z(Dim.rmc_1D.z+1:end,1:end)'*Lambda_1D.z(Dim.rmc_1D.z+1:end,1:end) , kron( I_m(Dim.rmc_1D.y) , I_m(m1) ) ) +...
         kron( I_m(m3) , kron( I_m(Dim.rmc_1D.y) , Lambda_1D.x(Dim.rmc_1D.x+1:end,1:end)'*Lambda_1D.x(Dim.rmc_1D.x+1:end,1:end) ) ) ,...   
         kron( I_m(Dim.rmc_1D.z) , kron( Lambda_1D.y(Dim.rmc_1D.y+1:end,1:end)'*Lambda_1D.y(Dim.rmc_1D.y+1:end,1:end) , I_m(m1) ) ) +...
         kron( I_m(Dim.rmc_1D.z) , kron( I_m(m2) , Lambda_1D.x(Dim.rmc_1D.x+1:end,1:end)'*Lambda_1D.x(Dim.rmc_1D.x+1:end,1:end) ) )  ...
         );
      
    % Corollary 2.
    % w = rand(3,1); w = w/norm(w);
    w = (1:3)'; w = w/norm(w);
    Lambda.Lambda_q = Lambda_3D.x'*Lambda_3D.x + Lambda_3D.y'*Lambda_3D.y + Lambda_3D.z'*Lambda_3D.z;
    Lambda_w = w(1)*Lambda_3D.x + w(2)*Lambda_3D.y + w(3)*Lambda_3D.z;
    
    norm_Pi_0     = sqrt(diag(Lambda.Lambda_q));
    Lambda.sqrt_Lambda_q = spdiags(norm_Pi_0,0,Dim.Nq,Dim.Nq);
    Lambda.inv_norm_Pi_0 = spdiags(1./norm_Pi_0, 0, Dim.Nq, Dim.Nq);
    norm_Pi_1     = sqrt(diag(Lambda.Lambda_q - Lambda_w'*Lambda_w));
    inv_norm_Pi_1 = spdiags(1./norm_Pi_1, 0, Dim.Nq, Dim.Nq);
    
    Pi.Pi_0  = [Lambda_3D.x;Lambda_3D.y;Lambda_3D.z]*Lambda.inv_norm_Pi_0;
    Pi.Pi_1  = -Lambda.r'*kron(w,I_m(Dim.Nq))*inv_norm_Pi_1;
    Pi.Pi_2  =  Lambda.r'*conj(Pi.Pi_1)*Lambda.inv_norm_Pi_0;
    Pi.D_Pi_1 = full([diag(Pi.Pi_1(1:Dim.Nq,:));diag(Pi.Pi_1(Dim.Nq+1:2*Dim.Nq,:));diag(Pi.Pi_1(2*Dim.Nq+1:3*Dim.Nq,:))]);
    Pi.D_Pi_2 = full([diag(Pi.Pi_2(1:Dim.Nq,:));diag(Pi.Pi_2(Dim.Nq+1:2*Dim.Nq,:));diag(Pi.Pi_2(2*Dim.Nq+1:3*Dim.Nq,:))]);
    
    Pi.Qq  = [       Pi.Pi_0,      Pi.Pi_1 ,      Pi.Pi_2  ];
    Pi.Pq  = [ conj(Pi.Pi_0),-conj(Pi.Pi_2), conj(Pi.Pi_1) ];
    Pi.Qqr = [      Pi.Pi_1 ,      Pi.Pi_2  ];
    Pi.Pqr = [-conj(Pi.Pi_2), conj(Pi.Pi_1) ];
    
    Lambda.Sigma.q  = blkdiag(sparse(Dim.Nq,Dim.Nq),Lambda.sqrt_Lambda_q,Lambda.sqrt_Lambda_q);
    Lambda.Sigma.qr = blkdiag(Lambda.sqrt_Lambda_q,Lambda.sqrt_Lambda_q);
    Lambda.Sigma.tilde_q = sqrt(Lambda.tilde_q);
    if Dim.Ntil == 0
        Lambda.Sigma.tilde_q_inv = [];
    else
        Lambda.Sigma.tilde_q_inv = spdiags(1./diag(Lambda.Sigma.tilde_q),0,Dim.Ntil,Dim.Ntil);
    end
    Lambda.Sigma.r = full([diag(Lambda.Sigma.qr);diag(Lambda.Sigma.tilde_q)]);
    Lambda.Sigma.inv_r = 1./Lambda.Sigma.r;    
    Lambda.Sigma.inv_r_sqr = Lambda.Sigma.inv_r.^2;
    O_x = DistOp.O_x(row, col);
    O_y = DistOp.O_y(row, col);
    O_z = DistOp.O_z(row, col);
    
    Pi.tilde = ...
    [ O_x               , Lambda_3D.Ctilde31, Lambda_3D.Ctilde21 ;
      Lambda_3D.Ctilde32, O_y               , Lambda_3D.Ctilde12 ;
      Lambda_3D.Ctilde23, Lambda_3D.Ctilde13, O_z               ];
    
    Pi.P  = [ I.hat, Pi.tilde*I.tilde*Lambda.Sigma.tilde_q_inv, I.Psi*Pi.Pq, Lambda.orig*I.tilde*Lambda.Sigma.tilde_q_inv ];
    Pi.Q  = [ I.Phi*Pi.Qq, I.tilde ];
    Pi.Pr = [ I.Psi*Pi.Pqr, Lambda.orig*I.tilde*Lambda.Sigma.tilde_q_inv ];
    Pi.Qr = [ I.Phi*Pi.Qqr, I.tilde ];
    Pi.Q0 =  I.Phi*Pi.Pi_0 ;
end