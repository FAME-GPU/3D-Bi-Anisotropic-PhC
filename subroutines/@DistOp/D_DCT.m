function vec = D_DCT(m)
    vec = [sqrt(1/2)*exp(1i*DistOp.theta(2*m, 1    )) ; 
                     exp(1i*DistOp.theta(2*m,(2:m)')) ];
    
end