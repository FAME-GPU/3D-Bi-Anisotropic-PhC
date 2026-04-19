function mtx = I_hat_3(row,col)
    mtx = kron( speye(col.z), kron( DistOp.eye_left(row.y,col.y), DistOp.eye_left(row.x,col.x)) );
end