function psf = MultiAperturePSF(xy_coords,aper_coords)
    % returns the point-spread function for the multi-aperture system
    n_apertures = size(aper_coords,1);
    n = 0;
    m = 0;
    v = 1;
    U = ones(1,n_apertures)/n_apertures; % psf is sum of all local modes
    
    psf = Basis_MixedAperture(xy_coords,n,m,v,U,aper_coords);
    
end