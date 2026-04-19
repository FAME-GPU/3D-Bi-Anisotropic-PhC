function N_int = IntCovarFormN_vertex(A, P, C)
% C = [C1,C2,C3] N传C   顶点插值 
for i = 1:3
    for j = 1:3
        if i ~= j
            tmp{i,j} = 1/4 * (2*P + C{i}) * A{i,j} * (2*P + C{j}');
        else
            tmp{i,j} = 1/2 * ( A{i,j} + (P + C{i}) * A{i,j} * (P + C{j})' );
        end
    end
end    
N_int = [tmp{1,1} ,tmp{1,2} ,tmp{1,3};
         tmp{2,1},tmp{2,2} ,tmp{2,3};
         tmp{3,1},tmp{3,2},tmp{3,3}];
end