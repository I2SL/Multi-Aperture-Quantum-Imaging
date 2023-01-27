function p = ModalProb_DirectImaging(s_x,s_y,X,Y,a_kx,a_ky)
    % Defines the probability distribution over the image plane for
    % given source positions s_x, s_y and aperture positions a_kx, a_ky.
    % In direct detection the measurement modes are delta functions 
    % on the image plane.
    
    dx = X(1,2) - X(1,1);
    dy = Y(2,1) - Y(1,1);
    p = zeros(numel(s_x),numel(X));
        
    % point spread function
    psf = @(x,y) MultiAperturePSF([x,y],[a_kx,a_ky]);
    
    % Shift psf for each source
    for k = 1:numel(s_x)
        p(k,:) = abs(psf(X(:)-s_x(k),Y(:)-s_y(k))).^2 *dx*dy;    
    end
    
end