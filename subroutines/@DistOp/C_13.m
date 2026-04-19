function mtx = C_13(row,col,delta,BC)
    mtx = kron( DistOp.K(row.z,delta.z,BC.z), kron( speye(row.y), speye(col.x) ) );
end