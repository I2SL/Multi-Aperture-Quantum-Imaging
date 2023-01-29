function GS_basis_pos_xy = Basis_GramSchmidt2(xy_coords,Kx,Ky,d2k,GS_basis_mom)
    % Computes the Gram-Schmidt basis in position space at the sample locations x,y
    % using an explicit Fourier Transform.
    
    x = xy_coords(:,1);
    y = xy_coords(:,2);
    
    
    % Perform manual FT in batches to avoid breaking memory limits for
    % large x,y vectors.
    n_modes = size(GS_basis_mom,2);
    n_pts = numel(x);
    batch_size = 200; % can be changed to accomodate memory limits of machine
    n_batches = ceil(n_pts/batch_size);
    remainder = rem(n_pts,batch_size);
    
    b_start = 1:batch_size:n_pts;
    b_end = batch_size:batch_size:n_pts;
    b_end = [b_end, b_end(end)+remainder];
    
    GS_basis_pos_xy = zeros(n_pts,n_modes);
    
    for k = 1:n_batches
        b = b_start(k):b_end(k); % batch indices
        
        % evaluation of FT at the location (x,y).
        FT_exp_xy = exp(1i * ( x(b).*Kx.' + y(b).*Ky.') );
        GS_basis_pos_xy(b,:) = 1/(2*pi) * d2k * FT_exp_xy * GS_basis_mom;
        
    end
        

end
