function mtx = DST_mtx(m)
    % m-by-m DST matrix 
    mtx = sqrt(2/(m+1)) * sin(    kron( (1:1:(  m))', DistOp.theta(m+1,2:m+1  ) ));
end