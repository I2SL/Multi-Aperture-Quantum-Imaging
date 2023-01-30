function dr_z = dr_FTZernike(r,theta,n,m)
    % derivative of the Fourier Transformed Zernike Polynomials
    % with respect to the radial parameter
    n_modes = numel(n);
    n_pts = numel(r);
    dr_z = ((-1).^(n/2)) ./(4*sqrt(pi)*sqrt(n+1)) .* FTzAngle(theta,m) .* ...
        (besselj(repmat(n-1,[n_pts,1]),repmat(r,[1,n_modes])) - besselj(repmat(n+3,[n_pts,1]),repmat(r,[1,n_modes])));
end
