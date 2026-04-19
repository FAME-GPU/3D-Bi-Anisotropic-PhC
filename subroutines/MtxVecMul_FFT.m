function vec_out = MtxVecMul_FFT(vec_in, Dim, istrans, S, D, P_idx_1, P_idx_2, P_idx_3, Type_1, Type_2, Type_3)
    
    global time_tmp
    % zss write
    if ~isfield(time_tmp, 'D')
        time_tmp.D = 0; time_tmp.S = 0; time_tmp.cat = 0;
        time_tmp.rs1 = 0; time_tmp.rs2 = 0; time_tmp.rs3 = 0;
        time_tmp.P1 = 0; time_tmp.P2 = 0; time_tmp.P3 = 0;
        time_tmp.fft1 = 0; time_tmp.fft2 = 0; time_tmp.fft3 = 0;
        time_tmp.ifft1 = 0; time_tmp.ifft2 = 0; time_tmp.ifft3 = 0;
        time_tmp.R1 = 0; time_tmp.R2 = 0; time_tmp.R3 = 0;
    end
    % zss end
    SIZE_1 = Dim.SIZE(1,:); SIZE_2 = Dim.SIZE(2,:); SIZE_3 = Dim.SIZE(3,:);
    N1 = length(P_idx_1); N2 = length(P_idx_2); N3 = length(P_idx_3);
    norm_c = sqrt(Dim.FFT_SIZE.x * Dim.FFT_SIZE.y * Dim.FFT_SIZE.z);
    
    if strcmp(istrans,'notranspose')
        %% MPV of D
tic;    x = D.*vec_in;                                    time_tmp.D = time_tmp.D + toc;
        % set 3D matrices x1->Y1, x2->Y2, x3->Y3
tic;    Y1 = reshape(x( (      1):(N1      ) ), SIZE_1);  time_tmp.rs1 = time_tmp.rs1 + toc;
tic;    Y2 = reshape(x( (   N1+1):(N1+N2   ) ), SIZE_2);  time_tmp.rs2 = time_tmp.rs2 + toc;
tic;    Y3 = reshape(x( (N1+N2+1):(N1+N2+N3) ), SIZE_3);  time_tmp.rs3 = time_tmp.rs3 + toc;
        %% MPV of P
tic;    Y1 = MtxVecMul_P(Y1, Type_1);                     time_tmp.P1 = time_tmp.P1 + toc;
tic;    Y2 = MtxVecMul_P(Y2, Type_2);                     time_tmp.P2 = time_tmp.P2 + toc;
tic;    Y3 = MtxVecMul_P(Y3, Type_3);                     time_tmp.P3 = time_tmp.P3 + toc;
        %%  IFFT step
tic;    Y1 = ifftn( Y1 )*norm_c;                          time_tmp.ifft1 = time_tmp.ifft1 + toc;
tic;    Y2 = ifftn( Y2 )*norm_c;                          time_tmp.ifft2 = time_tmp.ifft2 + toc;
tic;    Y3 = ifftn( Y3 )*norm_c;                          time_tmp.ifft3 = time_tmp.ifft3 + toc;
        %% MVP of R
tic;    Y1 = MtxVecMul_R(Y1, Type_1);  x1 = Y1(P_idx_1);  time_tmp.R1 = time_tmp.R1 + toc;
tic;    Y2 = MtxVecMul_R(Y2, Type_2);  x2 = Y2(P_idx_2);  time_tmp.R2 = time_tmp.R2 + toc;
tic;    Y3 = MtxVecMul_R(Y3, Type_3);  x3 = Y3(P_idx_3);  time_tmp.R3 = time_tmp.R3 + toc;
        %% MVP of S
tic;    x = cat(1,x1,x2,x3);                              time_tmp.cat = time_tmp.cat + toc;
tic;    vec_out = S.*x;                                   time_tmp.S = time_tmp.S + toc;
    elseif strcmp(istrans,'transpose')
        %% MVP of S'
tic;    x = conj(S).*vec_in;                              time_tmp.S = time_tmp.S + toc;
        % set 3D matrices x1->Y1, x2->Y2, x3->Y3
tic;    Y1 = reshape(x( (      1):(N1      ) ), SIZE_1);  time_tmp.rs1 = time_tmp.rs1 + toc;
tic;    Y2 = reshape(x( (   N1+1):(N1+N2   ) ), SIZE_2);  time_tmp.rs2 = time_tmp.rs2 + toc;
tic;    Y3 = reshape(x( (N1+N2+1):(N1+N2+N3) ), SIZE_3);  time_tmp.rs3 = time_tmp.rs3 + toc;
        %% MVP of R' 
tic;    Y1 = MtxVecMul_Rs(Y1, Type_1);                    time_tmp.R1 = time_tmp.R1 + toc;
tic;    Y2 = MtxVecMul_Rs(Y2, Type_2);                    time_tmp.R2 = time_tmp.R2 + toc;
tic;    Y3 = MtxVecMul_Rs(Y3, Type_3);                    time_tmp.R3 = time_tmp.R3 + toc;
        %%  FFT step
tic;    Y1 = fftn( Y1 )/norm_c;                           time_tmp.fft1 = time_tmp.fft1 + toc;
tic;    Y2 = fftn( Y2 )/norm_c;                           time_tmp.fft2 = time_tmp.fft2 + toc;
tic;    Y3 = fftn( Y3 )/norm_c;                           time_tmp.fft3 = time_tmp.fft3 + toc;
        %% MPV of P'
tic;    x1 = Y1(P_idx_1);                                 time_tmp.P1 = time_tmp.P1 + toc;
tic;    x2 = Y2(P_idx_2);                                 time_tmp.P2 = time_tmp.P2 + toc;
tic;    x3 = Y3(P_idx_3);                                 time_tmp.P3 = time_tmp.P3 + toc;
        %% MPV of D'
tic;    x = cat(1,x1,x2,x3);                              time_tmp.cat = time_tmp.cat + toc;
tic;    vec_out = conj(D).*x;                             time_tmp.D = time_tmp.D + toc;
    end
    clear x1 x2 x3
    
end

function X = MtxVecMul_P(X, Type)
    switch Type{1}
        case 'DCT'
            X = cat(1, X,...
                       zeros( size(X,1), size(X,2), size(X,3)));
        case 'DST'
            X = cat(1, zeros( 1, size(X,2), size(X,3)), X,...
                       zeros( size(X,1)+1, size(X,2), size(X,3)) );
    end
    switch Type{2}
        case 'DCT'
            X = cat(2, X,...
                       zeros( size(X,1), size(X,2), size(X,3)));
        case 'DST'
            X = cat(2, zeros( size(X,1), 1, size(X,3)), X,...
                       zeros( size(X,1), size(X,2)+1, size(X,3)));
    end
    switch Type{3}
        case 'DCT'
            X = cat(3, X,...
                       zeros( size(X,1), size(X,2), size(X,3)));
        case 'DST'
            X = cat(3, zeros( size(X,1), size(X,2), 1), X,...
                       zeros( size(X,1), size(X,2), size(X,3)+1));
    end
end


function Y = MtxVecMul_R(Y, Type)
    switch Type{1}
        case 'DCT'
            Y = Y + flip(Y,1);
        case 'DST'
            Z = flip(Y,1);
            Y = Y - circshift(Z,1,1);
    end
    switch Type{2}
        case 'DCT'
            Y = Y + flip(Y,2);
        case 'DST'
            Z = flip(Y,2);
            Y = Y - circshift(Z,1,2);
    end
    switch Type{3}
        case 'DCT'
            Y = Y + flip(Y,3);
        case 'DST'
            Z = flip(Y,3);
            Y = Y - circshift(Z,1,3);
    end
    clear Z
end

function X = MtxVecMul_Rs(X, Type)
    switch Type{1}
        case 'DCT'
            Z = flip(X,1);
            X = cat(1, X,...
                       Z);
        case 'DST'
            Z = flip(-X,1);
            X = cat(1, zeros( 1, size(X,2), size(X,3)), X,...
                       zeros( 1, size(X,2), size(X,3)), Z );
    end
    switch Type{2}
        case 'DCT'
            Z = flip(X,2);
            X = cat(2, X,...
                       Z);
        case 'DST'
            Z = flip(-X,2);
            X = cat(2, zeros( size(X,1), 1, size(X,3)), X,...
                       zeros( size(X,1), 1, size(X,3)), Z );
    end
    switch Type{3}
        case 'DCT'
            Z = flip(X,3);
            X = cat(3, X,...
                       Z);
        case 'DST'
            Z = flip(-X,3);
            X = cat(3, zeros( size(X,1), size(X,2), 1), X,...
                       zeros( size(X,1), size(X,2), 1), Z );
    end
    clear Z
end