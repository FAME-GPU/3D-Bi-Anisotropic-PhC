function mtx = K(m,delta,BC_info)
    switch BC_info.str
        case {'Periodic','periodic'}  
            mtx = DistOp.K_PBC(m,delta);
        case {'Quasi-Periodic','QP','quasi-periodic','quasiperiodic'}  
            mtx = DistOp.K_QPBC(m,BC_info.theta,delta);
        case {'Dirichlet','dirichlet','dirich','PEC'}
            mtx = DistOp.K_DBC(m,delta);
    end
end