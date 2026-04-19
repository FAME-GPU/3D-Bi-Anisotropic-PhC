function mtx = I_tilde_1(row,col)
    mtx = kron( speye(col.z), kron( speye(col.y), DistOp.eye_left(row.x,col.x)) );
end