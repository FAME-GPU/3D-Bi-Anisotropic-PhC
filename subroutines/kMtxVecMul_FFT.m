function vec_out = kMtxVecMul_FFT(vec_in, Dim, istrans, S, D, P_idx_1, P_idx_2, P_idx_3, Type_1, Type_2, Type_3)
    
    global x
    
    SIZE_1 = Dim.SIZE(1,:); SIZE_2 = Dim.SIZE(2,:); SIZE_3 = Dim.SIZE(3,:);
    N1 = length(P_idx_1); N2 = length(P_idx_2); N3 = length(P_idx_3);
    norm_c = sqrt(Dim.FFT_SIZE.x * Dim.FFT_SIZE.y * Dim.FFT_SIZE.z);
    
    if strcmp(istrans,'notranspose')
        %% MPV of D
        x = D.*vec_in;
        % set 3D matrices x1->Y1, x2->Y2, x3->Y3
        Y1 = reshape(x( (      1):(N1      ) ), SIZE_1);
        Y2 = reshape(x( (   N1+1):(N1+N2   ) ), SIZE_2);
        Y3 = reshape(x( (N1+N2+1):(N1+N2+N3) ), SIZE_3);
        %% MPV of P
%         tmp = MtxVecMul_P(Y1, Type_1);
        Y1 = kernel_MtxVecMul_P(Y1, Type_1);
        Y2 = MtxVecMul_P(Y2, Type_2);
        Y3 = MtxVecMul_P(Y3, Type_3);
        %%  IFFT step
        Y1 = ifftn( Y1 )*norm_c;
        Y2 = ifftn( Y2 )*norm_c;
        Y3 = ifftn( Y3 )*norm_c;
        %% MVP of R
        Y1 = MtxVecMul_R(Y1, Type_1); 
        Y2 = MtxVecMul_R(Y2, Type_2); 
        Y3 = MtxVecMul_R(Y3, Type_3); 
        x1 = Y1(P_idx_1);
        x2 = Y2(P_idx_2);
        x3 = Y3(P_idx_3);
        %% MVP of S
        x = cat(1,x1,x2,x3);
        vec_out = S.*x;
    elseif strcmp(istrans,'transpose')
        %% MVP of S'
        x = conj(S).*vec_in;
        % set 3D matrices x1->Y1, x2->Y2, x3->Y3
        Y1 = reshape(x( (      1):(N1      ) ), SIZE_1);
        Y2 = reshape(x( (   N1+1):(N1+N2   ) ), SIZE_2);
        Y3 = reshape(x( (N1+N2+1):(N1+N2+N3) ), SIZE_3);
        %% MVP of R'
        Y1 = MtxVecMul_Rs(Y1, Type_1);
        Y2 = MtxVecMul_Rs(Y2, Type_2);
        Y3 = MtxVecMul_Rs(Y3, Type_3);
        %%  FFT step
        Y1 = fftn( Y1 )/norm_c;
        Y2 = fftn( Y2 )/norm_c;
        Y3 = fftn( Y3 )/norm_c;
        %% MPV of P'
        x1 = Y1(P_idx_1);
        x2 = Y2(P_idx_2);
        x3 = Y3(P_idx_3);
        %% MPV of D'
        x = cat(1,x1,x2,x3);
        vec_out = conj(D).*x;
    end
    clear x1 x2 x3
    
end

function X = MtxVecMul_P(X, Type)
    switch Type{1}
        case 'DCT'
            X = cat(1, X,...
                       zeros( size(X,1), size(X,2), size(X,3), 'gpuArray'));
        case 'DST'
            X = cat(1, zeros( 1, size(X,2), size(X,3), 'gpuArray'), X,...
                       zeros( size(X,1)+1, size(X,2), size(X,3), 'gpuArray') );
    end
    switch Type{2}
        case 'DCT'
            X = cat(2, X,...
                       zeros( size(X,1), size(X,2), size(X,3), 'gpuArray'));
        case 'DST'
            X = cat(2, zeros( size(X,1), 1, size(X,3), 'gpuArray'), X,...
                       zeros( size(X,1), size(X,2)+1, size(X,3), 'gpuArray'));
    end
    switch Type{3}
        case 'DCT'
            X = cat(3, X,...
                       zeros( size(X,1), size(X,2), size(X,3), 'gpuArray'));
        case 'DST'
            X = cat(3, zeros( size(X,1), size(X,2), 1, 'gpuArray'), X,...
                       zeros( size(X,1), size(X,2), size(X,3)+1, 'gpuArray'));
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
            X = cat(1, zeros( 1, size(X,2), size(X,3), 'gpuArray'), X,...
                       zeros( 1, size(X,2), size(X,3), 'gpuArray'), Z );
    end
    switch Type{2}
        case 'DCT'
            Z = flip(X,2);
            X = cat(2, X,...
                       Z);
        case 'DST'
            Z = flip(-X,2);
            X = cat(2, zeros( size(X,1), 1, size(X,3), 'gpuArray'), X,...
                       zeros( size(X,1), 1, size(X,3), 'gpuArray'), Z );
    end
    switch Type{3}
        case 'DCT'
            Z = flip(X,3);
            X = cat(3, X,...
                       Z);
        case 'DST'
            Z = flip(-X,3);
            X = cat(3, zeros( size(X,1), size(X,2), 1, 'gpuArray'), X,...
                       zeros( size(X,1), size(X,2), 1, 'gpuArray'), Z );
    end
    clear Z
end

function X = kernel_MtxVecMul_P(Y,Type)
    switch Type{1}
        case 'DFT'
            mode(1) = 0;
            size_out(1) = size(Y,1);
        case 'DCT'
            mode(1) = 1;
            size_out(1) = 2*size(Y,1);
        case 'DST'
            mode(1) = 2;
            size_out(1) = 2*(size(Y,1)+1);
    end
    switch Type{2}
        case 'DFT'
            mode(2) = 0;
            size_out(2) = size(Y,2);
        case 'DCT'
            mode(2) = 1;
            size_out(2) = 2*size(Y,2);
        case 'DST'
            mode(2) = 2;
            size_out(2) = 2*(size(Y,2)+1);
    end
    switch Type{3}
        case 'DFT'
            mode(3) = 0;
            size_out(3) = size(Y,3);
        case 'DCT'
            mode(3) = 1;
            size_out(3) = 2*size(Y,3);
        case 'DST'
            mode(3) = 2;
            size_out(3) = 2*(size(Y,3)+1);
    end

    k = parallel.gpu.CUDAKernel('MtxVecMul_P.ptx','MtxVecMul_P.cu');
    k.GridSize = [size(Y,2),size(Y,3)];
    k.ThreadBlockSize = size(Y,1);
    X = zeros(size_out(1),size_out(2),size_out(3),'gpuArray') + 1i*zeros(size_out(1),size_out(2),size_out(3),'gpuArray');
%     X = complex(parallel.gpu.GPUArray.zeros(size_out,'double')); 
    X = feval(k,X,Y,size(Y,1),size(Y,2),size(Y,3),mode(1),mode(2),mode(3));
    
end