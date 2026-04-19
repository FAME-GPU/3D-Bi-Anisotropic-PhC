function mtx = C_b1(row,col,delta,BC)
    mtx = kron( speye(row.z), kron( speye(row.y), DistOp.K(row.x,delta.x,BC.x) ) );
end