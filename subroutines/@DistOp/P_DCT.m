function mtx = P_DCT(m)
    mtx = [speye(m);sparse(m,m)];
end