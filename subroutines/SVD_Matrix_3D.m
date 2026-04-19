function [Lambda_3D, R, R_idx1, R_idx2, P, P_idx, D, S, Type] = SVD_Matrix_3D(Lambda_1D, V, U, Dim)
    I_m       = @(m)   speye(m);
    I_mn      = @(m,n) [sparse(m-n,n);I_m(n)];
    
    n1 = Dim.n1; n2 = Dim.n2; n3 = Dim.n3;
    m1 = Dim.m1; m2 = Dim.m2; m3 = Dim.m3;
    
    Lambda_3D.C13 = kron( Lambda_1D.z , kron(     I_m(n2) ,     I_m(m1) ) );
    Lambda_3D.C12 = kron(     I_m(n3) , kron( Lambda_1D.y ,     I_m(m1) ) );
    Lambda_3D.C23 = kron( Lambda_1D.z , kron(     I_m(m2) ,     I_m(n1) ) );
    Lambda_3D.C21 = kron(     I_m(n3) , kron(     I_m(m2) , Lambda_1D.x ) );
    Lambda_3D.C32 = kron(     I_m(m3) , kron( Lambda_1D.y ,     I_m(n1) ) );
    Lambda_3D.C31 = kron(     I_m(m3) , kron(     I_m(n2) , Lambda_1D.x ) );
  
    Lambda_3D.x = kron( I_m(m3) , kron( I_m(m2) , Lambda_1D.x(Dim.rmc_1D.x+1:end,1:end) ) ) ;
    Lambda_3D.y = kron( I_m(m3) , kron( Lambda_1D.y(Dim.rmc_1D.y+1:end,1:end), I_m(Dim.m1)  ) ) ;   
    Lambda_3D.z = kron( Lambda_1D.z(Dim.rmc_1D.z+1:end,1:end), kron( I_m(Dim.m2) , I_m(Dim.m1)  ) ) ; 
    
    Lambda_3D.Ctilde31 = kron( I_mn(n3,m3) , kron( I_m(n2) , Lambda_1D.x(Dim.rmc_1D.x+1:end,1:end)' ) ) ;
    Lambda_3D.Ctilde21 = kron( I_m(n3) , kron( I_mn(n2,m2) , Lambda_1D.x(Dim.rmc_1D.x+1:end,1:end)' ) ) ;
    Lambda_3D.Ctilde32 = kron( I_mn(n3,m3) , kron( Lambda_1D.y(Dim.rmc_1D.y+1:end,1:end)' , I_m(n1) ) ) ;
    Lambda_3D.Ctilde12 = kron( I_m(n3) , kron( Lambda_1D.y(Dim.rmc_1D.y+1:end,1:end)' , I_mn(n1,m1) ) ) ;
    Lambda_3D.Ctilde23 = kron( Lambda_1D.z(Dim.rmc_1D.z+1:end,1:end)' , kron( I_mn(n2,m2) , I_m(n1) ) ) ;
    Lambda_3D.Ctilde13 = kron( Lambda_1D.z(Dim.rmc_1D.z+1:end,1:end)' , kron( I_m(n2) , I_mn(n1,m1) ) ) ;

    
    
    R.Psi_1 = kron( V.z.R, kron( V.y.R, U.x.R) );
    R.Psi_2 = kron( V.z.R, kron( U.y.R, V.x.R) );
    R.Psi_3 = kron( U.z.R, kron( V.y.R, V.x.R) );
    R.Psi = blkdiag(R.Psi_1, R.Psi_2, R.Psi_3);
    P.Psi_1 = kron( V.z.P, kron( V.y.P, U.x.P) );
    P.Psi_2 = kron( V.z.P, kron( U.y.P, V.x.P) );
    P.Psi_3 = kron( U.z.P, kron( V.y.P, V.x.P) );
    P.Psi = blkdiag(P.Psi_1, P.Psi_2, P.Psi_3);
    D.Psi_1 = kron( V.z.D, kron( V.y.D, U.x.D));
    D.Psi_2 = kron( V.z.D, kron( U.y.D, V.x.D));
    D.Psi_3 = kron( U.z.D, kron( V.y.D, V.x.D));
    D.Psi = [D.Psi_1; D.Psi_2; D.Psi_3];
    
    R.start_idx.Psi = [ U.x.start_idx, V.y.start_idx, V.z.start_idx;
                        V.x.start_idx, U.y.start_idx, V.z.start_idx;
                        V.x.start_idx, V.y.start_idx, U.z.start_idx ];
    R.start_idx.Phi = [ V.x.start_idx, U.y.start_idx, U.z.start_idx;
                        U.x.start_idx, V.y.start_idx, U.z.start_idx;
                        U.x.start_idx, U.y.start_idx, V.z.start_idx ];
    
    R.Phi_1 = kron( U.z.R, kron( U.y.R, V.x.R) );
    R.Phi_2 = kron( U.z.R, kron( V.y.R, U.x.R) );
    R.Phi_3 = kron( V.z.R, kron( U.y.R, U.x.R) );
    R.Phi = blkdiag(R.Phi_1, R.Phi_2, R.Phi_3);
    P.Phi_1 = kron( U.z.P, kron( U.y.P, V.x.P) );
    P.Phi_2 = kron( U.z.P, kron( V.y.P, U.x.P) );
    P.Phi_3 = kron( V.z.P, kron( U.y.P, U.x.P) );
    P.Phi = blkdiag(P.Phi_1, P.Phi_2, P.Phi_3);
    D.Phi_1 = kron( U.z.D, kron( U.y.D, V.x.D));
    D.Phi_2 = kron( U.z.D, kron( V.y.D, U.x.D));
    D.Phi_3 = kron( V.z.D, kron( U.y.D, U.x.D));
    D.Phi = [D.Phi_1; D.Phi_2; D.Phi_3];
    
    P.start_idx.Psi = [ U.x.start_idx, V.y.start_idx, V.z.start_idx;
                        V.x.start_idx, U.y.start_idx, V.z.start_idx;
                        V.x.start_idx, V.y.start_idx, U.z.start_idx ];
    P.start_idx.Phi = [ V.x.start_idx, U.y.start_idx, U.z.start_idx;
                        U.x.start_idx, V.y.start_idx, U.z.start_idx;
                        U.x.start_idx, U.y.start_idx, V.z.start_idx ];
                 
    Dim.FFT_SIZE.x = size(V.x.R,2);
    Dim.FFT_SIZE.y = size(V.y.R,2);
    Dim.FFT_SIZE.z = size(V.z.R,2);
    
    Type.Phi_1 = {V.x.type, U.y.type, U.z.type};
    Type.Phi_2 = {U.x.type, V.y.type, U.z.type};
    Type.Phi_3 = {U.x.type, U.y.type, V.z.type};
    Type.Psi_1 = {U.x.type, V.y.type, V.z.type};
    Type.Psi_2 = {V.x.type, U.y.type, V.z.type};
    Type.Psi_3 = {V.x.type, V.y.type, U.z.type};
    
    S.Psi_1 = kron(V.z.S, kron(V.y.S, U.x.S));
    S.Psi_2 = kron(V.z.S, kron(U.y.S, V.x.S));
    S.Psi_3 = kron(U.z.S, kron(V.y.S, V.x.S));
    S.Phi_1 = kron(U.z.S, kron(U.y.S, V.x.S));
    S.Phi_2 = kron(U.z.S, kron(V.y.S, U.x.S));
    S.Phi_3 = kron(V.z.S, kron(U.y.S, U.x.S));
    S.Psi = [S.Psi_1; S.Psi_2; S.Psi_3];
    S.Phi = [S.Phi_1; S.Phi_2; S.Phi_3];
    
    R_idx111.Psi_1 = find( kron( V.z.R_idx1, kron( V.y.R_idx1, U.x.R_idx1) ) == 1 );
    R_idx112.Psi_1 = find( kron( V.z.R_idx1, kron( V.y.R_idx1, U.x.R_idx2) ) == 1 );
    R_idx121.Psi_1 = find( kron( V.z.R_idx1, kron( V.y.R_idx2, U.x.R_idx1) ) == 1 );
    R_idx122.Psi_1 = find( kron( V.z.R_idx1, kron( V.y.R_idx2, U.x.R_idx2) ) == 1 );
    R_idx211.Psi_1 = find( kron( V.z.R_idx2, kron( V.y.R_idx1, U.x.R_idx1) ) == 1 );
    R_idx212.Psi_1 = find( kron( V.z.R_idx2, kron( V.y.R_idx1, U.x.R_idx2) ) == 1 );
    R_idx221.Psi_1 = find( kron( V.z.R_idx2, kron( V.y.R_idx2, U.x.R_idx1) ) == 1 );
    R_idx222.Psi_1 = find( kron( V.z.R_idx2, kron( V.y.R_idx2, U.x.R_idx2) ) == 1 );
    R_idx111.Psi_2 = find( kron( V.z.R_idx1, kron( U.y.R_idx1, V.x.R_idx1) ) == 1 );
    R_idx112.Psi_2 = find( kron( V.z.R_idx1, kron( U.y.R_idx1, V.x.R_idx2) ) == 1 );
    R_idx121.Psi_2 = find( kron( V.z.R_idx1, kron( U.y.R_idx2, V.x.R_idx1) ) == 1 );
    R_idx122.Psi_2 = find( kron( V.z.R_idx1, kron( U.y.R_idx2, V.x.R_idx2) ) == 1 );
    R_idx211.Psi_2 = find( kron( V.z.R_idx2, kron( U.y.R_idx1, V.x.R_idx1) ) == 1 );
    R_idx212.Psi_2 = find( kron( V.z.R_idx2, kron( U.y.R_idx1, V.x.R_idx2) ) == 1 );
    R_idx221.Psi_2 = find( kron( V.z.R_idx2, kron( U.y.R_idx2, V.x.R_idx1) ) == 1 );
    R_idx222.Psi_2 = find( kron( V.z.R_idx2, kron( U.y.R_idx2, V.x.R_idx2) ) == 1 );
    R_idx111.Psi_3 = find( kron( U.z.R_idx1, kron( V.y.R_idx1, V.x.R_idx1) ) == 1 );
    R_idx112.Psi_3 = find( kron( U.z.R_idx1, kron( V.y.R_idx1, V.x.R_idx2) ) == 1 );
    R_idx121.Psi_3 = find( kron( U.z.R_idx1, kron( V.y.R_idx2, V.x.R_idx1) ) == 1 );
    R_idx122.Psi_3 = find( kron( U.z.R_idx1, kron( V.y.R_idx2, V.x.R_idx2) ) == 1 );
    R_idx211.Psi_3 = find( kron( U.z.R_idx2, kron( V.y.R_idx1, V.x.R_idx1) ) == 1 );
    R_idx212.Psi_3 = find( kron( U.z.R_idx2, kron( V.y.R_idx1, V.x.R_idx2) ) == 1 );
    R_idx221.Psi_3 = find( kron( U.z.R_idx2, kron( V.y.R_idx2, V.x.R_idx1) ) == 1 );
    R_idx222.Psi_3 = find( kron( U.z.R_idx2, kron( V.y.R_idx2, V.x.R_idx2) ) == 1 );
    
    P_idx.Psi_1 = find( kron( V.z.P_idx, kron( V.y.P_idx, U.x.P_idx) ) == 1 );
    P_idx.Psi_2 = find( kron( V.z.P_idx, kron( U.y.P_idx, V.x.P_idx) ) == 1 );
    P_idx.Psi_3 = find( kron( U.z.P_idx, kron( V.y.P_idx, V.x.P_idx) ) == 1 );
    
    R_idx1.Phi_1 = find( kron( U.z.R_idx1, kron( U.y.R_idx1, V.x.R_idx1) ) == 1 );
    R_idx2.Phi_1 = find( kron( U.z.R_idx2, kron( U.y.R_idx2, V.x.R_idx2) ) == 1 );
    R_idx1.Phi_2 = find( kron( U.z.R_idx1, kron( V.y.R_idx1, U.x.R_idx1) ) == 1 );
    R_idx2.Phi_2 = find( kron( U.z.R_idx2, kron( V.y.R_idx2, U.x.R_idx2) ) == 1 );
    R_idx1.Phi_3 = find( kron( V.z.R_idx1, kron( U.y.R_idx1, U.x.R_idx1) ) == 1 );
    R_idx2.Phi_3 = find( kron( V.z.R_idx2, kron( U.y.R_idx2, U.x.R_idx2) ) == 1 );
    P_idx.Phi_1 = find( kron( U.z.P_idx, kron( U.y.P_idx, V.x.P_idx) ) == 1 );
    P_idx.Phi_2 = find( kron( U.z.P_idx, kron( V.y.P_idx, U.x.P_idx) ) == 1 );
    P_idx.Phi_3 = find( kron( V.z.P_idx, kron( U.y.P_idx, U.x.P_idx) ) == 1 );
end