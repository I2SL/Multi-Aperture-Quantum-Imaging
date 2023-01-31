function aper_coords  = Golay9(R)
    % Generates the aperture coordinates for a Golay-9 telescope array.
    % R is the radius of the outer most circle on which the apertures
    % reside.
    R1 = R * 1/3 * ones(3,1); 
    R2 = R * 2/3 * ones(3,1);
    R3 = R * 3/3 * ones(3,1);
    rho = [R1; R2; R3];   % golay-9 radii 

    tri = linspace(0,(2*pi)*2/3,3)' + pi/2;
    a1 = (2*pi/3)*0/3 + tri;
    a2 = (2*pi/3)*1/3 + tri;
    a3 = (2*pi/3)*2/3 + tri;
    phi = [a1; a2; a3];   % golay-9 angles 
    
    [a_kx, a_ky] = pol2cart(phi,rho);
    aper_coords = [a_kx,a_ky];
end
