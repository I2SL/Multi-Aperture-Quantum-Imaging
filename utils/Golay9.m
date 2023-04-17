function aper_coords  = Golay9(R)
    % Generates the aperture coordinates for a Golay-9 telescope array.
    % R is the radius of the outer most circle on which 3 of the apertures
    % reside. The inner circles on which the remaining 6 apertures reside
    % have radii 2/3 R and 1/3 R.
    
    R1 = R * 1/3 * ones(3,1);   % radial coordinates for inner circle of apertures
    R2 = R * 2/3 * ones(3,1);   % radial coordinates for middle circle of apertures
    R3 = R * 3/3 * ones(3,1);   % radial coordinates for outer circle of apertures
    rho = [R1; R2; R3];   % golay-9 radii 

    tri = linspace(0,(2*pi)*2/3,3)' + pi/2;
    a3 = (2*pi/3)*0/3 + tri;    % angular coordinates for inner circle
    a2 = (2*pi/3)*1/3 + tri;    % angular coordiantes for middle circle
    a1 = (2*pi/3)*2/3 + tri;    % angular coordinates for outer circle
    phi = [a1; a2; a3];   % golay-9 angles 
    
    [a_kx, a_ky] = pol2cart(phi,rho);
    aper_coords = [a_kx,a_ky];
end
