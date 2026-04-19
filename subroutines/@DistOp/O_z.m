function mtx = O_z(row,col)
    mtx = sparse(row.x*row.y*col.z, col.x*col.y*row.z);
end