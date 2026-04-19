function y = fncB(x, fnc, mtx_B)
%UNTITLED2 此处显示有关此函数的摘要
y1 = [fnc.Pr(x(1:end/2)); fnc.Qr(x(end/2+1:end))];
y2 = [mtx_B.M_int * y1(1:end/2) + mtx_B.J_int' * y1(end/2+1:end); ...
      mtx_B.J_int * y1(1:end/2) + mtx_B.N_int * y1(end/2+1:end)];
y =  [fnc.Prs(y2(1:end/2)); fnc.Qrs(y2(end/2+1:end))];
end