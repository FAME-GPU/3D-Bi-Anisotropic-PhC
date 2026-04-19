function mtx = K_QPBC(m,theta,delta)
    % forward difference matrix with QPBC
    mtx = sparse([1:m,1:m]', [1:m,2:m,1]', [-ones(m,1); ones(m-1,1); exp(1i*theta)], m, m) / delta;   
end
