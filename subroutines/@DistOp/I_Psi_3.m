function mtx = I_Psi_3(row,col)
    mtx = kron( speye(col.z), kron( DistOp.eye_right(row.y,col.y), DistOp.eye_right(row.x,col.x)) );
end