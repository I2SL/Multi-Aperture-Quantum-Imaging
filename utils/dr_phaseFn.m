function dr_B = dr_phaseFn(xy_coords,v,U,aper_coords) 
    % partial derivative of the phase function with respect to the radial
    % coordinate
    
    % phase introduced by each shifted sub-aperture
    phase = exp( 1i * aper_coords * xy_coords');
    
    % direction projections
    xy_hat_coords = xy_coords./vecnorm(xy_coords,2,2);
    aper_cos = 1i * aper_coords * xy_hat_coords';
 
    % derivative of modal interference term with respect to source separation
    dr_B = (U(v,:) * (aper_cos .* phase)).';
end

