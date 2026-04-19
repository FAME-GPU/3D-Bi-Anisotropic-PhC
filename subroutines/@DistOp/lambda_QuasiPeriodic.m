function vec = lambda_QuasiPeriodic(m,theta,delta)
    % lambda_QuasiPeriodic : singular values of the m-by-m matrix Kqp(m,theta)
    vec = (exp(1i*(2*DistOp.theta(m,(1:m)')+theta/m)) - 1)/delta; 
end

        