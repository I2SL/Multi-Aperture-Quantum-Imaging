function psi_nm_v = Basis_MixedAperture(xy_coords,n,m,v,U,aper_coords)
    % Returns the the mixed aperture-local modal basis. The chosen aperture-local
    % modes for this study are the fourier transform of the Zernike
    % polynomials - a basis of orthonormal functions defined over the
    % unit disk.
    

    % xy_coords : Nx2 matrix with cartesian coordinate pairs at which to
    % evaluate the basis functions
    % n,m : 1xM vectors containing Zernike mode indices
    % v   : 1xM vectors containing mixed mode indices 
    % U   : KxK unitary mixing matrix defining how the local modes will be mixed
    % aper_coords: Kx2 matrix with cartesian coordinate pairs for the aperture positions
    %--------------
    % psi_nm - [NxM]

    
    
    % phase from multi-aperture mixing
    B = phaseFn(xy_coords,v,U,aper_coords);
    
    % mixed mode
    [theta,r] = cart2pol(xy_coords(:,1),xy_coords(:,2));
    psi_nm_v = B .* FTZernike(r,theta,n,m);
end