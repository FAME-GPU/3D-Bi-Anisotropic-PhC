function mtx = C_b3(row,col,delta,BC)
    mtx = kron( DistOp.K(row.z,delta.z,BC.z), kron( speye(row.y), speye(row.x) ) );
end