function [ e1, e2, e3, h1, h2, h3, SFD, omega ] = Eigen_Restoration( EV, lambda, Dim, mtx_B, fnc_Qr, fnc_Pr, fnc_tildeM, Sigma_r, NFSEP_type, opt )
    
    M = mtx_B.M_int;
    N = mtx_B.N_int;

    %% Restore eigenpairs
    omega = zeros(size(EV,2),1);
    for j = 1 : size(EV,2)
        if NFSEP_type == 0
            x = EV(:,j) / norm(EV(:,j));
            w = 1./sqrt(lambda(j));
            er = Sigma_r.\x;
            hr = -1i*fnc_tildeM(Sigma_r.*er);
        elseif NFSEP_type ==1
            x = EV(      1:end/2,j) / norm(EV(      1:end/2,j));
            y = EV(end/2+1:end  ,j) / norm(EV(end/2+1:end  ,j));
            w = 1i./lambda(j);
            er = Sigma_r.\y;
            hr = Sigma_r.\x;
        end
        omega(j) = w;
        h  = -1i*(M*fnc_Pr(Sigma_r.*er))/w;
        e  =  1i*(N*fnc_Qr(Sigma_r.*hr))/w;
        
        e1(:,j) = e(                1:Dim.Ne1        );
        e2(:,j) = e(        Dim.Ne1+1:Dim.Ne1+Dim.Ne2);
        e3(:,j) = e(Dim.Ne1+Dim.Ne2+1:end            );
        h1(:,j) = h(                1:Dim.Nh1        );
        h2(:,j) = h(        Dim.Nh1+1:Dim.Nh1+Dim.Nh2);
        h3(:,j) = h(Dim.Nh1+Dim.Nh2+1:end            );
        
        Dim = opt.comp.Dim;
        n1 = Dim.n1; m1 = Dim.m1;
        n2 = Dim.n2; m2 = Dim.m2;
        n3 = Dim.n3; m3 = Dim.m3;  

        row.x = n1; row.y = n2; row.z = n3;
        col.x = m1; col.y = m2; col.z = m3;
        delta.x = 1; delta.y = 1; delta.z = 1;
        
        BC = opt.phys.BC;
        
        R{1} = DistOp.I_R1(row,col); R{2} = DistOp.I_R2(row,col); R{3} = DistOp.I_R3(row,col);
        S{1} = DistOp.I_S1(row,col); S{2} = DistOp.I_S2(row,col); S{3} = DistOp.I_S3(row,col);
        Cd{1} = -DistOp.C_d1(row,col,delta,BC); Cd{2} = -DistOp.C_d2(row,col,delta,BC); Cd{3} = -DistOp.C_d3(row,col,delta,BC);
        Cb{1} =  DistOp.C_b1(row,col,delta,BC); Cb{2} =  DistOp.C_b2(row,col,delta,BC); Cb{3} =  DistOp.C_b3(row,col,delta,BC);
        
        e1_tmp = (Cd{1} + 2*R{1})*e1(:,j)/2;
        e2_tmp = (Cd{2} + 2*R{2})*e2(:,j)/2;
        e3_tmp = (Cd{3} + 2*R{3})*e3(:,j)/2;
        SFD(:,j) =  abs(e1_tmp.*e1_tmp + e2_tmp.*e2_tmp + e3_tmp.*e3_tmp) ;
    end

    %% Reordering
    [omega,idx] = sort(omega);
    e1 = e1(:,idx); e2 = e2(:,idx); e3 = e3(:,idx);
    h1 = h1(:,idx); h2 = h2(:,idx); h3 = h3(:,idx);
    SFD = SFD(:,idx);
end