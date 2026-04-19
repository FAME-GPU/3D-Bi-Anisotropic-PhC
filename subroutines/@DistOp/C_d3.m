function mtx = C_d3(row,col,delta,BC)
    mtx = kron( -DistOp.K(row.z,delta.z,BC.z)', kron( speye(col.y), speye(col.x) ) );
end