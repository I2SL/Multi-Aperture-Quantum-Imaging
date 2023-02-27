function [Kx,Ky,d2k] = ApertureKxKy(aper_coords,aper_rads,samp_dense)
% aper_coords - [Mx2] matrix of (kx, ky) coordinates for the each aperture centroid
% aper_radius -  [Mx1] vector of aperture radii 
% samp_dense  - [double] number samples per unit length in k-space

% returns the support of the aperture in K-space coordinates
num_apertures = size(aper_coords,1);

% differential element
dkx = 1/samp_dense; dky = dkx;
d2k = dkx*dky;

% K-space coordinates inside each circular aperture with uniform sampling density
Kx = [];
Ky = [];

for j = 1:num_apertures
    % bounding box for each sub-aperture
    r_j = aper_rads(j);             % radius of j'th aperture
    s_j = ceil(r_j*samp_dense);     % number of samples along diameter of aperture
    s_j = s_j + ~mod(s_j,2);        % make sure the number of samples is odd
    [kx,ky] = meshgrid(linspace(-r_j,r_j,s_j)); 
    
    % centered sub-aperture coordinates
    in_j = kx.^2 + ky.^2 < r_j^2;
    
    kx_j = kx(in_j);
    ky_j = ky(in_j);
    
    Kx = [Kx; kx_j + aper_coords(j,1)];
    Ky = [Ky; ky_j + aper_coords(j,2)];

end

end





%{
function [Kx,Ky,d2k] = ApertureKxKy(aper_coords,subap_sampling)
% returns the support of a system of circular apertures of diameter 1 in
% k_space coordinates.
num_apertures = size(aper_coords,1);

% differential element
dkx = 1/subap_sampling; dky = dkx;
d2k = dkx*dky;

% K-space coordinates inside each circular aperture with uniform sampling
% density
Kx = [];
Ky = [];


for j = 1:num_apertures
    % bounding box for each sub-aperture
    [kx,ky] = meshgrid(linspace(-1,1,subap_sampling)); 
        
    % centered sub-aperture coordinates
    in_j = kx.^2 + ky.^2 < (.5)^2;
    
    kx_j = kx(in_j);
    ky_j = ky(in_j);
    
    Kx = [Kx; kx_j + aper_coords(j,1)];
    Ky = [Ky; ky_j + aper_coords(j,2)];

end



end
%}

