function psi_nm_v = Basis_MixedAperture(xy_coords,n,m,v,U,aperture)
    % Returns the the mixed aperture-local modal basis. The chosen aperture-local
    % modes for this study are the fourier transform of the Zernike
    % polynomials - a basis of orthonormal functions defined over the
    % unit disk.
    

    % xy_coords : Nx2 matrix with cartesian coordinate pairs at which to
    % evaluate the basis functions
    % n,m : 1xM vectors containing Zernike mode indices
    % v   : 1xM vectors containing mixed mode indices 
    % U   : KxK unitary mixing matrix defining how the local modes will be mixed
    % aperture : Kx3 matrix with cartesian coordinate pairs for the
    %           aperture positions and the radii
    %--------------
    % psi_nm - [NxM]
       

    [theta,r] = cart2pol(xy_coords(:,1),xy_coords(:,2));
    
    % get unique mixing indices
    vv = unique(v);
    % get unique zernike indices
    z = unique([n;m]','rows')';
    nn = z(1,:); mm = z(2,:); 
    
    
    % mixed mode
    N = size(xy_coords,1); % number of coordinate query points
    Z = size(z,2);         % number of zernike modes
    K = size(aperture,1);
    psi_nm_k = zeros(N,Z,K);

    % get the local aperture modes with phase coefficient introduced by
    % shift position of each sub-apertures
    for k = 1:K
        h_k = aperture(k,1:2); % aperture position
        r_k = aperture(k,3);   % aperture radius        
        psi_nm_k(:,:,k) = r_k*exp(1i*xy_coords*h_k.').*FTZernike(r_k*r,theta,nn,mm);
    end
    
    % apply mixing transform
    psi_nm_v = pagemtimes(U,permute(psi_nm_k,[3,2,1]));
    psi_nm_v = permute(psi_nm_v,[3,2,1]);
    psi_nm_v = psi_nm_v(:,:,vv);
    psi_nm_v = reshape(psi_nm_v,[N,Z*numel(vv)]);   
end

%{
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
    % psi_nm - [NxMxK]

    % phase from multi-aperture mixing
    B = phaseFn(xy_coords,v,U,aper_coords);
    
    % mixed mode
    [theta,r] = cart2pol(xy_coords(:,1),xy_coords(:,2));
    psi_nm_v = B .* FTZernike(r,theta,n,m);
      
end
%}
