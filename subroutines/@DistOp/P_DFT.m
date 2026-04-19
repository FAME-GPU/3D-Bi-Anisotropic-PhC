function mtx = P_DFT(m)
    mtx = sparse((2:m)',(1:m-1)',ones(m-1,1),m,m-1);
end