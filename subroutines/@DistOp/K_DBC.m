function mtx = K_DBC(m,delta)
    % forward difference matrix with DBC
    mtx = sparse([2:m,1:m-1]', [1:m-1,1:m-1]', [-ones(m-1,1); ones(m-1,1)], m, m-1) / delta;
end

        