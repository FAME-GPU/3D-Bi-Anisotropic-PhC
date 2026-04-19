function mtx = K_PBC(m,delta)
    % forward difference matrix with PBC
    mtx = sparse([1:m,1:m]', [1:m,2:m,1]', [-ones(m,1); ones(m-1,1); 1], m, m) / delta; 
end