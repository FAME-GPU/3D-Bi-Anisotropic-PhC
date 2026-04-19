function mtx = R_DST(m)
    mtx = [sparse(m-1,1),speye(m-1),sparse(m-1,1),-DistOp.eye_flip(m-1)];
end