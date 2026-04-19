function vec = MtxVecMul_Br(x, fnc1, fnc2, mtx )
    vec = fnc1(x);
    vec = mtx.*vec;
    vec = fnc2(vec);
end