function mtx = O_x(row,col)
    mtx = sparse(col.x*row.y*row.z, row.x*col.y*col.z);
end