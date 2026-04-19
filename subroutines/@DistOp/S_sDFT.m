function vec = S_sDFT(m,theta) 
    vec = exp(1i*(0:1:(m-1))'*theta/m);
end