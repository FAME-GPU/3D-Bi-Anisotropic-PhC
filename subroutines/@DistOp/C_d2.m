function mtx = C_d2(row,col,delta,BC)
    mtx = kron( speye(col.z), kron( -DistOp.K(row.y,delta.y,BC.y)', speye(col.x) ) );
end