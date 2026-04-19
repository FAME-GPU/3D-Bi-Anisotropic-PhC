function vec = MtxVecMul_Ar(x, mtx, fnc )
    vec = mtx.*x;
    vec = fnc(vec);
    vec = mtx.*vec;
end