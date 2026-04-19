function mtx = C_d1(row,col,delta,BC)
    mtx = kron( speye(col.z), kron( speye(col.y), -DistOp.K(row.x,delta.x,BC.x)' ) );
end