function f = ctsIFT_2D(x,y,kx,ky,d2k,F_tilde) 
    % Computes the Continuous 2D Inverse Fourier Transform of the function 
    % F_tilde defined at kx,ky at the locations x,y.
    %
    %
    % x,y       - [N,1] vectors denoting position-space coordinates to evaluate
    %           the 2D IFT at.
    %
    % kx,ky     - [M,1] vectors denoting the momentum-space coordinates over
    %           which F_tilde is defined.
    %
    % d2k       - scalar denoting the 2d differential element for integration
    %           over momentum space.
    %
    % F_tilde   - [M,K] matrix defining the set of K functions to perform the inverse
    %           transform on. Each column of F_tilde is a unique function.
    %           F_tilde(i,j) is the j'th function evaluted at momentum-space coordinate (kx(i),ky(i)).
    % --------------------------------------------------------------------
    % f         - [N,K] matrix containing the set of K IFT'd functions.
    %           Each column of f is a unique function.
    %           f(i,j) is the j'th IFT'd function evaluated at
    %           position-space coordinate (x(i),y(i));
    
    

    
    % generate computation batches to handle memory limits
    n_pts = numel(x);
    batch_size = 200;
    n_batches = ceil(n_pts/batch_size);
    remainder = rem(n_pts,batch_size);
    
    % batch starting and stopping indices
    b_start = 1:batch_size:n_pts;
    b_end = [batch_size:batch_size:n_pts,batch_size*(n_batches-1)+remainder];
    
    num_fs = size(F_tilde,2);
    
    f = zeros(n_pts,num_fs);
    
    % exploit GPU
    if gpuDeviceCount("available") > 0
        f = gpuArray(f);
        kx = gpuArray(kx);
        ky = gpuArray(ky);
        x = gpuArray(x);
        y = gpuArray(y);
    end

    for i = 1:n_batches
        b = b_start(i):b_end(i); % batch indices

        % manual evaluation of FT at the location (x,y).
        FT_exp_xy = exp(1i * ( x(b).*kx.' + y(b).*ky.') );
        f(b,:) = FT_exp_xy * F_tilde;
        
    end
    
    f = gather(f) * 1/(2*pi) * d2k;

end

