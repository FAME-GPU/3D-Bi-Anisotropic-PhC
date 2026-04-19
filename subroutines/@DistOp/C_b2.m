function mtx = C_b2(row,col,delta,BC)
    mtx = kron( speye(row.z), kron( DistOp.K(row.y,delta.y,BC.y), speye(row.x) ) );
end