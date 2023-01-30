function [Kx,Ky,d2k] = ApertureKxKy(aper_coords,circle_sampling)
% returns the support of the aperture in K-space coordinates
num_apertures = size(aper_coords,1);

% bounding box for a sub-aperture
[kx,ky] = meshgrid(linspace(-1,1,circle_sampling)); 

% differential element
dkx = (kx(1,2)-kx(1,1));
dky = (ky(2,1)-ky(1,1));
d2k = dkx * dky;

% centered sub-aperture support coordinates
in_circle = kx.^2 + ky.^2 < 1;
kx_circle = kx(in_circle)'; 
ky_circle = ky(in_circle)';
num_cicor = nnz(in_circle);      % number of kx,ky coordinates pairs inside each sub-aperture


% generate k-space coordinates that correspond to the support of the
% multi-aperture system
Kx = zeros(num_cicor*num_apertures,1);
Ky = zeros(num_cicor*num_apertures,1);
for j = 1:num_apertures
    j_rng = (num_cicor*(j-1)+1):(num_cicor*j);
    Kx(j_rng) = kx_circle + aper_coords(j,1);
    Ky(j_rng) = ky_circle + aper_coords(j,2); 
end

end