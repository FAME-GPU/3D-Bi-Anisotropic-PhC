function vec = MtxVecMul_Qr(x, istrans, fnc, mtx )
    if strcmp(istrans,'notranspose')
        vec = fnc( mtx*x );
    elseif strcmp(istrans,'transpose')
        vec = mtx*fnc( x );
    end
end