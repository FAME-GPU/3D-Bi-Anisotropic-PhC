function [vec1, vec2, sgn] = Ridx_DCT(m)
    vec1 = [ones(m,1);zeros(m,1)];
    vec2 = [zeros(m,1);ones(m,1)];
    sgn = 1;