function mtx = I_Psi_1(row,col)
    mtx = kron( DistOp.eye_right(row.z,col.z), kron( DistOp.eye_right(row.y,col.y), speye(col.x)) );
end