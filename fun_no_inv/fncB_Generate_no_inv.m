function fncB = fncB_Generate_no_inv(fnc,mtx_B,y)
% 先求逆再插值的特征值问题的左端矩阵生成
% x1 = [fnc.Qr(y(end/2+1:end));fnc.Pr(y(1:end/2))];
% x2 = [mtx_B.N_int * x1(1:end/2) + mtx_B.J_int * x1(end/2+1:end) ;mtx_B.J_int' * x1(1:end/2) + mtx_B.M_int * x1(end/2+1:end)];
% fncB = [fnc.Prs(x2(end/2+1:end));fnc.Qrs(x2(1:end/2))];
x1 = [fnc.Pr(y(1:end/2));fnc.Qr(y(end/2+1:end))];
x2 = [mtx_B.M_int * x1(1:end/2) + mtx_B.J_int' * x1(end/2+1:end) ;mtx_B.J_int * x1(1:end/2) + mtx_B.N_int * x1(end/2+1:end)];
fncB = [fnc.Prs(x2(1:end/2));fnc.Qrs(x2(end/2+1:end))];
end