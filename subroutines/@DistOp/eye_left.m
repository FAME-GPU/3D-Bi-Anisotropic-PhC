function mtx = eye_left(m,n)
    mtx = [speye(m-n);sparse(n,m-n)];
end