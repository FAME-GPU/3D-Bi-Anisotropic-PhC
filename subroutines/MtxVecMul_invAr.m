function vec = MtxVecMul_invAr(x, mtx, fnc )
    vec = mtx.\x;
    vec = fnc(vec);
    vec = mtx.\vec;
end