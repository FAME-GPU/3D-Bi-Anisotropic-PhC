function mtx = C_12(row,col,delta,BC)
    mtx = kron( speye(row.z), kron( DistOp.K(row.y,delta.y,BC.y), speye(col.x) ) );
end