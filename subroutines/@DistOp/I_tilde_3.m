function mtx = I_tilde_3(row,col)
    mtx = kron( DistOp.eye_left(row.z,col.z), kron( speye(col.y), speye(col.x)) );
end