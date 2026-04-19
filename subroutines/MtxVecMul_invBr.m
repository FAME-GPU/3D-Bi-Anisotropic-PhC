function vec = MtxVecMul_invBr(x, fnc, linearsolver )
    if isempty(linearsolver)~=1
        x0    = linearsolver.x0_fnc(length(x));
        switch linearsolver.solver
            case 'pcg'
            vec = pcg(fnc, x, linearsolver.tol, linearsolver.itmax, [], [], x0);
            case 'bicgstabl'
            vec = bicgstabl(fnc, x, linearsolver.tol, linearsolver.itmax, [], [], x0);
        end
    end
end