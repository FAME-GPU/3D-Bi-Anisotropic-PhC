function mtx = I_Phi_3(row,col)
    mtx = kron( DistOp.eye_right(row.z,col.z), kron( speye(col.y), speye(col.x)) );
end