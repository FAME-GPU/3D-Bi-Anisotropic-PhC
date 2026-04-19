function mtx = C_23(row,col,delta,BC)
    mtx = kron( DistOp.K(row.z,delta.z,BC.z), kron( speye(col.y), speye(row.x) ) );
end