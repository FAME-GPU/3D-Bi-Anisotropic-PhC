function mtx = C_21(row,col,delta,BC)
    mtx = kron( speye(row.z), kron( speye(col.y), DistOp.K(row.x,delta.x,BC.x) ) );
end