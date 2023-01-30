function [Kx,Ky,d2k,GS_basis_mom] = genGramSchmidtBasis_mom(max_order,aper_coords,circle_sampling)
%function [Kx,Ky,GS_basis_mom] = genGramSchmidtBasis2(n_max,aper_radii,circle_sampling)
% Generates the 2D PSF-Adapted Gram-Schmidt basis for any aperture function
% Kx,Ky - coordinate pairs in K-space over the aperture support (in units of sub-aperture radius)
% A - The normalized aperture function defined over Kx Ky
% n_max - max polynomial order

num_modes = (max_order+1)^2;
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


% instantiate discretized GS basis function matrix and polynomial
% coefficient chart
num_supcrd = numel(Kx);                         % number of support coordinates
GS_basis_mom = zeros(num_supcrd,num_modes);     % discretized GS mode matrix
%poly_coeff = zeros(n_max+1,n_max+1,num_modes);  % array of P_nm(ikx ,iky) expansion polynomial coefficients for each basis function

mode = 1;
for n = 0:max_order
    for m = 0:max_order
    
        % moments for candidate polynomial
        ikxn = (1i*Kx).^n;
        ikym = (1i*Ky).^m;
        
        % candidate polynomial
        pol_nm = ikxn.*ikym;
        
        % projections onto existing basis functions
        proj_coeff = GS_basis_mom'*pol_nm*d2k;
        
        % remove projections
        phi_nm = pol_nm - GS_basis_mom*proj_coeff;
        
        % normalize
        phi_nm = phi_nm / sqrt( phi_nm'*phi_nm * d2k);
        
        % add new basis element
        GS_basis_mom(:,mode) = phi_nm;
        
        % index mode
        mode = mode+1;
    end   
end


% visualize
%figure(1)
%GS_i = abs(GS_basis_mom(:,36)).^2;
%scatter3(Kx,Ky,GS_i,2,'filled')


end

