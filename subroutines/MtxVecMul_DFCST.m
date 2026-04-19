function x = MtxVecMul_DFCST(x, Dim_1, Dim_2, Dim_3, Type_1, Type_2, Type_3, type, fnc, S)
    N1 = prod(Dim_1);
    N2 = prod(Dim_2);
    N3 = prod(Dim_3);

    if strcmp(type,'notranspose')
tic;    Y1 = reshape(x( (      1):(N1      ) ), Dim_1);
tic;    Y2 = reshape(x( (   N1+1):(N1+N2   ) ), Dim_2);
tic;    Y3 = reshape(x( (N2+N1+1):(N1+N2+N3) ), Dim_3);
    
tic;    Y1 = IDFCST(Y1, Type_1, fnc);
tic;    Y2 = IDFCST(Y2, Type_2, fnc);
tic;    Y3 = IDFCST(Y3, Type_3, fnc);
        
tic;    x = S.*[Y1(:); Y2(:); Y3(:)];
    elseif strcmp(type,'transpose')
tic;    x = conj(S).*x;
        
tic;    Y1 = reshape(x( (      1):(N1      ) ), Dim_1);
tic;    Y2 = reshape(x( (   N1+1):(N1+N2   ) ), Dim_2);
tic;    Y3 = reshape(x( (N2+N1+1):(N1+N2+N3) ), Dim_3);
        
tic;    Y1 = DFCST(Y1, Type_1, fnc);
tic;    Y2 = DFCST(Y2, Type_2, fnc);
tic;    Y3 = DFCST(Y3, Type_3, fnc);
        
tic;    x = [Y1(:); Y2(:); Y3(:)];
    end
    
end

function Y = DFCST(Y, Type, fnc)
    Y = eval(['fnc.',Type{1},'_x(Y)']);
    Y = eval(['fnc.',Type{2},'_y(Y)']);
    Y = eval(['fnc.',Type{3},'_z(Y)']);
end
function Y = IDFCST(Y, Type, fnc)
    Y = eval(['fnc.I',Type{1},'_x(Y)']);
    Y = eval(['fnc.I',Type{2},'_y(Y)']);
    Y = eval(['fnc.I',Type{3},'_z(Y)']);
end
    