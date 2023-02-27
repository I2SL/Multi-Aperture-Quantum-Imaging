function GS_basis_mom = genGramSchmidtBasis_mom(max_order,Kx,Ky,d2k)
% Generates the 2D PSF-Adapted Gram-Schmidt basis for any binary aperture function
% max_order - the maximum order of GS bases desired (number of derivatives to perform in x and y)
% Kx,Ky - coordinate pairs in K-space over the aperture support (in units of sub-aperture radius)
% d2k - The

num_modes = (max_order+1)^2;


% instantiate discretized GS basis function matrix and polynomial
% coefficient chart
num_supcrd = numel(Kx);                         % number of support coordinates
GS_basis_mom = zeros(num_supcrd,num_modes);     % discretized GS mode matrix

if gpuDeviceCount('available')>0
    GS_basis_mom = gpuArray(GS_basis_mom);
    Kx = gpuArray(Kx);
    Ky = gpuArray(Ky);
end

tic
mode = 1;
for n = 0:max_order
    for m = 0:max_order
    
        % moments for candidate polynomial
        ikxn = (1i*Kx).^n;
        ikym = (1i*Ky).^m;
        
        % candidate polynomial
        poly_nm = ikxn.*ikym;
        
        % projections onto existing basis functions
        proj_coeff = GS_basis_mom'*poly_nm*d2k;
        
        % remove projections
        phi_nm = poly_nm - GS_basis_mom*proj_coeff;
        
        % normalize
        phi_nm = phi_nm / sqrt( phi_nm'*phi_nm * d2k);
        
        % add new basis element
        GS_basis_mom(:,mode) = phi_nm;
        
        % index mode
        mode = mode+1;
    end   
end

GS_basis_mom = gather(GS_basis_mom);
disp('GS modes created')
toc
end

