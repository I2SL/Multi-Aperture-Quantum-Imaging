function coords = PlusAperture(n,R)
    % returns the multi-aperture centroid coordinates for an aperture that
    % looks like a '+' sign. The sub-apertures are located on the x and y
    % axes only. The maximum radial extent of the aperture coordinates is D.
    % The number of apertures must be chosen to have symmetric
    % configuration. Therefore the number of apertures n must be either a
    % multiple of 4 or one greater than a multiple of 4 (for an aperture to
    % be located at the origin).

    assert(n >0)
    assert(rem(n,4) < 1e-15 || (rem(n,4)-1) < 1e-15)

    coords = [];
    if rem(n,4) == 1
        coords = [0,0];
    end
    
    m = 2*floor(n/4) + rem(n,4);
    
    arm = R*linspace(-1,1,m)';
    arm(arm == 0) = [];
    
    coords = [coords;
              [zeros(numel(arm),1),arm];
              [arm,zeros(numel(arm),1)]
              ];
end