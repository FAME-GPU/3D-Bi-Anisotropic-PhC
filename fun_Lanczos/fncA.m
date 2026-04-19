function y = fncA(x, mtx_A)
%UNTITLED2 此处显示有关此函数的摘要
y = [1i * mtx_A.Lambda_curl.Sigma.inv_r .* x(end/2+1:end);...
   - 1i * mtx_A.Lambda_curl.Sigma.inv_r .* x(1:end/2)];
end