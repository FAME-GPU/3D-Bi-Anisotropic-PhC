function mtx = R_DCT(m)
    mtx = [speye(m),DistOp.eye_flip(m)];
end