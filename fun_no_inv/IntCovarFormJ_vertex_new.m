function J_int = IntCovarFormJ_vertex_new(A, I, C)
% C = [C1,C2,C3] N传C   顶点插值 
for i = 1:3
    for j = 1:3
        k = setdiff([1 2 3], j);
         tmp{i,j} = 1/8 * (2*I + C{i}) * A{i,j} * (2*I + C{k(2)}') * (2*I + C{k(1)}');
    end
end    
J_int = [tmp{1,1},tmp{1,2},tmp{1,3};
         tmp{2,1},tmp{2,2},tmp{2,3};
         tmp{3,1},tmp{3,2},tmp{3,3}];
end