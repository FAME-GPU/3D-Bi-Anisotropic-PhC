function M_int = IntCovarFormM_vertex_new(A, I, C)
% C = [C1,C2,C3] M传C   顶点插值 
for i = 1:3 
    k = setdiff([1 2 3], i);
    for j = 1:3
        if i == j            
            tmp{i,i} = 1/4*(A{i,i} + (I + C{k(1)}) * A{i,i} * (I + C{k(1)}') ...
                       + (I + C{k(2)}) * A{i,i} * (I + C{k(2)}') ...
                       + (I + C{k(1)}) * (I + C{k(2)}) * A{i,i} * (I + C{k(2)}') * (I + C{k(1)}') );
        else
            l = setdiff([1 2 3], j);
            tmp{i,j} = 1/16 * (2*I + C{k(1)}) * (2*I + C{k(2)}) * A{i,j} * (2*I + C{l(2)}') * (2*I + C{l(1)}');
            % tmp{i,j} = 1/16 * (2*I + C{k(1)}) * (2*I + C{k(2)}) * A{i,j} * (2*I + C{k(2)}') * (2*I + C{k(1)}');
        end
    end
end    
M_int = [tmp{1,1} ,tmp{1,2} ,tmp{1,3};
         tmp{2,1},tmp{2,2} ,tmp{2,3};
         tmp{3,1},tmp{3,2},tmp{3,3}];
end