function mtx = I_hat_1(row,col)
    mtx = kron( DistOp.eye_left(row.z,col.z), kron( DistOp.eye_left(row.y,col.y), speye(col.x)) );
end