function p = ModalProb_DirectImaging(src_coords,X,Y,aperture)
    % Defines the probability distribution over the image plane for
    % given source positions s_x, s_y and aperture positions a_kx, a_ky.
    % In direct detection the measurement modes are delta functions 
    % on the image plane.
    s_x = src_coords(:,1);
    s_y = src_coords(:,2);
    
    dx = X(1,2) - X(1,1);
    dy = Y(2,1) - Y(1,1);
    d2r = dx*dy;
    p = zeros(numel(s_x),numel(X));
        
    % point spread function
    psf = @(x,y) MultiAperturePSF([x,y],aperture);
    
    % Shift psf for each source
    for k = 1:numel(s_x)
        p(k,:) = abs(psf(X(:)-s_x(k),Y(:)-s_y(k))).^2 *d2r;    
    end
    
end