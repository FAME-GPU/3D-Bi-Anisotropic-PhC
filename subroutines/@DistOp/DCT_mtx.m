function mtx = DCT_mtx(m)
    % m-by-m DCT matrix 
    mtx = [ sqrt(1/m) * cos(0.5*kron( (1:2:(2*m-1))', DistOp.theta(m,1    ) )),...  
            sqrt(2/m) * cos(0.5*kron( (1:2:(2*m-1))', DistOp.theta(m,2:m  ) )) ];
end