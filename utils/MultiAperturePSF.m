function psf = MultiAperturePSF(xy_coords,aperture)
    n = 0;
    m = 0;
    v = 1;
    
    ap_num = size(aperture,1);
    U = ones(1,ap_num)/sqrt(ap_num); % psf is sum of all local modes
    
    % PSF is sum of phase shifted local 0 modes with uniform mixing
    psf = Basis_MixedAperture(xy_coords,n,m,v,U,aperture); 
end

%{
function psf = MultiAperturePSF(xy_coords,aper_coords)
    
    % convert to polar coordinates
    [theta,r] = cart2pol(xy_coords(:,1),xy_coords(:,2));
    
    % returns the point-spread function for the multi-aperture system
    n_apertures = size(aper_coords,1);
    n = 0;
    m = 0;
    v = 1;
    U = ones(1,n_apertures)/sqrt(n_apertures); % psf is sum of all local modes
    
    
    % Mode and modulation
    Z = FTZernike(r,theta,n,m);                         % single aperture PSF
    B = phaseFn(xy_coords,v,U,aper_coords);             % phases introduced by multiple apertures
    
    psf = Z.*B;
end
%}