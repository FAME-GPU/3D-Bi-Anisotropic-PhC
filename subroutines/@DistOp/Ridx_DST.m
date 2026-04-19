function [vec1, vec2, sgn] = Ridx_DST(m)
    vec1 = [0;ones(m-1,1);zeros(m,1)];
    vec2 = [zeros(m,1);0;ones(m-1,1)];
    sgn = -1;