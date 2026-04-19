function mtx_B = Interpolation_Material_Matrix(opt, mtx_B)
    Dim = opt.comp.Dim;
    n1 = Dim.n1; m1 = Dim.m1;
    n2 = Dim.n2; m2 = Dim.m2;
    n3 = Dim.n3; m3 = Dim.m3;  

    %% Construct covariance form of parameter matrices
    row.x = n1; row.y = n2; row.z = n3;
    col.x = m1; col.y = m2; col.z = m3;
    delta.x = 1; delta.y = 1; delta.z = 1;
    
    BC = opt.phys.BC;
   
    R{1} = DistOp.I_R1(row,col); R{2} = DistOp.I_R2(row,col); R{3} = DistOp.I_R3(row,col);
    S{1} = DistOp.I_S1(row,col); S{2} = DistOp.I_S2(row,col); S{3} = DistOp.I_S3(row,col);
    Cd{1} = -DistOp.C_d1(row,col,delta,BC); Cd{2} = -DistOp.C_d2(row,col,delta,BC); Cd{3} = -DistOp.C_d3(row,col,delta,BC);
    Cb{1} =  DistOp.C_b1(row,col,delta,BC); Cb{2} =  DistOp.C_b2(row,col,delta,BC); Cb{3} =  DistOp.C_b3(row,col,delta,BC);
    
    
    % mtx_B.N_int = IntCovarForm(mtx_B.N, R, Cd);
    % mtx_B.M_int = IntCovarForm(mtx_B.M, S, Cb);
    
    %% 手性
    % N: \varepsilon; M: \mu;  L: \zeta;  J: \xi.
    mtx_B.N_int = Int_simple(mtx_B.N);
    mtx_B.M_int = Int_simple(mtx_B.M);
    mtx_B.L_int = Int_simple(mtx_B.L);
    mtx_B.J_int = Int_simple(mtx_B.J);
    % 辅助计算
    mtx_B.M_inv_int = Int_simple(mtx_B.M_inv);
    mtx_B.Phi_inv_int = Int_simple(mtx_B.Phi_inv);
end


function N_int = IntCovarForm(A, P, C)
    for i = 1:3
        for j = 1:3
            if i == j
                tmp{i,i} = (P{i}'*A{i,i}*P{i} + (P{i} + C{i})'*A{i,i}*(P{i} + C{i}))/2;
            else
                tmp{i,j} = (2*P{i} + C{i})'*A{i,j}*(2*P{j} + C{j})/4;
            end
        end
    end
    
    N_int = [tmp{1,1},tmp{1,2},tmp{1,3};
             tmp{2,1},tmp{2,2},tmp{2,3};
             tmp{3,1},tmp{3,2},tmp{3,3}];
end

function N_int = Int_simple(tmp)
    N_int = [tmp{1,1},tmp{1,2},tmp{1,3};
             tmp{2,1},tmp{2,2},tmp{2,3};
             tmp{3,1},tmp{3,2},tmp{3,3}];
end
