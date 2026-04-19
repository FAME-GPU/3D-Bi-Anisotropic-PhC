function mtx = C_31(row,col,delta,BC)
    mtx = kron( speye(col.z), kron( speye(row.y),DistOp.K(row.x,delta.x,BC.x) ) );
end