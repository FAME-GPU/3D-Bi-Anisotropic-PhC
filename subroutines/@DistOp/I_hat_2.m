function mtx = I_hat_2(row,col)
    mtx = kron( DistOp.eye_left(row.z,col.z), kron( speye(col.y), DistOp.eye_left(row.x,col.x)) );
end