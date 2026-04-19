function x = gDST( x, n, dim )

    
    if isempty(n)
            m = 1 + size(x,dim);
    else
        m = 1 + n;
    end
    norm_fact = sqrt(2/m);

    idx_x = 1;
    for i = 1:length(size(x))
        if i == dim
            tmp = [0;ones(m-1,1,'gpuArray');zeros(m,1,'gpuArray')];
        else
            tmp = ones(size(x,i),1,'gpuArray');
        end
        idx_x = kron(tmp,idx_x);
    end
    idx_valid = find(idx_x==1);
    
    size_y = size(x);
    size_y(dim) = 2*m;
    
    y = zeros(size_y,'gpuArray');
    y(idx_valid) = x;

    y = fft(y,2*m,dim);
    
    x = reshape(y(idx_valid),size(x,1),size(x,2),size(x,3));
    x = norm_fact*(-1)*imag(x);
end