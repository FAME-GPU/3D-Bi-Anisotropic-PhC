function y = fncinvB(x, B_op, tol, maxit)
% 返回fncB\x
% [y, flag, relres, iter] = pcg(@(z)fncB(z, fnc, mtx_B), x, tol, maxit);
[y, ~] = pcg(B_op, x, tol, maxit);
end