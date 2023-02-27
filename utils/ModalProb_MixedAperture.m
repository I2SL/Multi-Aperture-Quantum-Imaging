function p = ModalProb_MixedAperture(xy_coords,n,m,v,U,aperture)
    % returns a discrete PMF over the 2D modes which constitutes the photon modal detection
    % probabilities for photons emitted by a point source at the position x,y

    % correlation function
    correlation_fn = corrFn_MixedAperture(xy_coords,n,m,v,U,aperture);
    
    p = abs(correlation_fn).^2;    % probability is norm squared of correlation function
end

%{
function p = ModalProb_MixedAperture(xy_coords,n,m,v,U,aper_coords,A_tot)
    % returns a discrete PMF over the 2D modes which constitutes the photon modal detection
    % probabilities for photons emitted by a point source at the position x,y

    % correlation function
    correlation_fn = corrFn_MixedAperture(xy_coords,n,m,v,U,aper_coords,A_tot);
    
    p = abs(correlation_fn).^2;    % probability is norm squared of correlation function
end
%}