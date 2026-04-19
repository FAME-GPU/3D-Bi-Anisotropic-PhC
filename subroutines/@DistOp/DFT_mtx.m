function mtx = DFT_mtx(m)
    % F_m : m-by-m DFT matrix
    mtx = sqrt(1/m) * exp(2i*kron( (0:1:(  m-1))', DistOp.theta(m,1:m) ));
end