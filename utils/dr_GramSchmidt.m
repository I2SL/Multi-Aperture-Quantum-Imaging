function dr_Gamma_pos_xy = dr_GramSchmidt(xy_coords,Kx,Ky,d2k,GS_basis_mom)
    % Computes the derivative of the the Gram-Schmidt basis
    % with respect to radial coordinate using the differentiation property
    % of the Fourier Transform.
    
    % get the angular coordinate
    x = xy_coords(:,1);
    y = xy_coords(:,2);
    [theta,~] = cart2pol(x,y);

    % Perform manual FT in batches to avoid breaking memory limits for
    % large x,y vectors.
    n_modes = size(GS_basis_mom,2);
    n_pts = numel(x);
    batch_size = 200; % can be changed to accomodate memory limits of machine
    n_batches = ceil(n_pts/batch_size);
    remainder = rem(n_pts,batch_size);
    
    b_start = 1:batch_size:n_pts;
    b_end = [batch_size:batch_size:n_pts,batch_size*(n_batches-1)+remainder];
    
    
    dr_Gamma_pos_xy = zeros(n_pts, n_modes);
    
    
    % exploit GPU
    if gpuDeviceCount("available") > 0
        dr_gamma_pos_xy = gpuArray(dr_Gamma_pos_xy);
        kx = gpuArray(Kx);
        ky = gpuArray(Ky);
        x = gpuArray(x);
        y = gpuArray(y);
    end

    
    for  k = 1:n_batches
        b = b_start(k):b_end(k); % batch indices
        
        % evaluation of d/dr FT at the location (x,y)
        FT_exp_xy = exp(1i * ( x(b).*Kx.' + y(b).*Ky.') );
        dr_FT_exp_xy = 1i * ( cos(theta).*Kx.' + sin(theta).*Ky.').*FT_exp_xy;
        dr_Gamma_pos_xy(b,:) =  dr_FT_exp_xy * GS_basis_mom;
    end
    
    dr_Gamma_pos_xy = 2*gather(dr_Gamma_pos_xy) * 1/(2*pi) * d2k;  % not sure why factor of 2 is needed here but it makes the probability in the 0th mode = 1 for 0 source displacement as desired.
    
end