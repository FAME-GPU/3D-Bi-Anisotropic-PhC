function vec = lambda_Dirichlet(m,delta)
    % lambda_Dirichlet : singular values of the m-by-(m-1) matrix Kd(m)
    vec = 2*sin(0.5*DistOp.theta(m,(2:m)')) / delta; 
end