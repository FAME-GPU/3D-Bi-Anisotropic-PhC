function mtx = O_y(row,col)
    mtx = sparse(row.x*col.y*row.z, col.x*row.y*col.z);
end