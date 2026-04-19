function mtx_B = Interpolation_Material_Matrix_Mult_Point_vertex(opt, mtx_B, type)
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
    %Cd=[Cd{1},Cd{2},Cd{3}]=[C1^*,C2^*,C3^*]
    %Cb=[Cb{1},Cb{2},Cb{3}]=[C1,C2,C3]
    
    if strcmp(type, 'bianisotripic') 
        mtx_B.N_int = IntCovarFormN_vertex(mtx_B.N, R{1}, Cb);
        mtx_B.N_int = 1/2*( mtx_B.N_int + mtx_B.N_int'); 
        mtx_B.M_int = IntCovarFormM_vertex_new(mtx_B.M, R{1}, Cb); 
        mtx_B.M_int = 1/2*( mtx_B.M_int + mtx_B.M_int');   
        mtx_B.J_int = IntCovarFormJ_vertex_new(mtx_B.J, R{1}, Cb);
    end
end
