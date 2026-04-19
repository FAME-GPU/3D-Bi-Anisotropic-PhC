function vec = lambda_Periodic(m,delta)
    % lambda_Periodic : singular values of the m-by-(m matrix Kp(m)
    vec = DistOp.lambda_QuasiPeriodic(m,0,delta);
end

        