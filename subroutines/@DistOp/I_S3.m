function mtx = I_S3(row,col)
    mtx = kron( DistOp.eye_right(row.z,col.z), kron( speye(row.y), speye(row.x) ) );
end