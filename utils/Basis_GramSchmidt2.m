function GS_basis_pos_xy = Basis_GramSchmidt2(x,y,Kx,Ky,d2k,GS_basis_mom)
    % Manually evaluates the GS basis in position space at the sample locations x,y
    % as this is a fourier transform over a non-uniform point set.
    
    
    % Perform manual FT in batches to avoid breaking memory limits
    n_modes = size(GS_basis_mom,2);
    n_pts = numel(x);
    batch_size = 200;
    n_batches = ceil(n_pts/batch_size);
    remainder = rem(n_pts,batch_size);
    
    b_start = 1:batch_size:n_pts;
    b_end = batch_size:batch_size:n_pts;
    b_end = [b_end, b_end(end)+remainder];
    
    GS_basis_pos_xy = zeros(n_pts,n_modes);
    
    for k = 1:n_batches
        b = b_start(k):b_end(k); % batch indices
        
        % manual evaluation of FT at the location (x,y).
        FT_exp_xy = exp(1i * ( x(b).*Kx.' + y(b).*Ky.') );
        GS_basis_pos_xy(b,:) = 1/(2*pi) * d2k * FT_exp_xy * GS_basis_mom;
        
    end
        
    
    %{
    pre-vectorization implementation
    
    
    %d2k = (Kx(1,2)-Kx(1,1))*(Ky(2,1)-Ky(1,1));
 
    % vectorize the modes
    GS_basis_mom = reshape(GS_basis_mom, [1,numel(Kx),size(GS_basis_mom,3)]);
    
    FT_exp_xy = exp(1i * ( x.*Kx(:)' + y.*Ky(:)' ) ); % exponential terms in fourier transform
    GS_basis_pos_xy = 1/(2*pi) * sum( GS_basis_mom .* FT_exp_xy ,2) .* d2k;

    GS_basis_pos_xy = squeeze(GS_basis_pos_xy);
    %}
end
