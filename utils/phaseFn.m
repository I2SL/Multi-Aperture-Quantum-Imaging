function B = phaseFn(xy_coords,v,U,aper_coords)
    % Mixed-aperture phase function that arises from interfering the modes
    
    % phase introduced by each shifted sub-aperture
    phase = exp( 1i * aper_coords * xy_coords');
    
    % aperture mixing
    B = (U(v,:) * phase).';
end
