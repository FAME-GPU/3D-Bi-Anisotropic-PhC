function mtx = I_Phi_1(row,col)
    mtx = kron( speye(col.z), kron( speye(col.y), DistOp.eye_right(row.x,col.x)) );
end