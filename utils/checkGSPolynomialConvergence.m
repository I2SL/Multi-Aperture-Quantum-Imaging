function checkGSPolynomialConvergence(aper_coords,n_max)
% check GS coefficient convergence
ap_dims = [51,101:100:1001]; % 51 is a pre-estimate

poly_coeff_prev = zeros(n_max+1,n_max+1,(n_max+1)^2); 
for i = 1:numel(ap_dims)
    ap_dim = ap_dims(i);

    % Gram-Schmidt Modes
    [aperture,Kx,Ky] = ApertureConfig(aper_coords(:,1),aper_coords(:,2),ap_dim);
    [poly_coeff,~] = GSPolynomials(Kx,Ky,aperture,n_max);

    % compute max difference
    eps = 1e-1; % threshold for distinguishing relevant coefficients from irrelevant ones
    idx = abs(poly_coeff)>eps;
    abs_error(i,:) = abs(poly_coeff(idx) - poly_coeff_prev(idx)); % absolute error (fractional error doesnt matter because small-value polynomials don't contribute to GS basis)
    
    poly_coeff_prev = poly_coeff;
end

figure
plot(ap_dims(2:end),abs_error(2:end,:))
title('Gram-Schmidt Polynomial Coefficient Convergence')
xlabel('Aperture Plane Samples (1D)')
ylabel('Running Absolute Difference')

end


function [poly_coeff,poly_basis] = GSPolynomials(Kx,Ky,A,n_max)
% Generates the 2D PSF-Adapted Gram-Schmidt momentum-space polynomials
% for any aperture function
% 
% Kx,Ky - a meshgrid of K-space coordinates
% A - The normalized aperture function defined over Kx Ky
% n_max - max polynomial order
% ap_dim - the square aperture plane 1D dimensionality
% ip_dim - the square image plane1D  dimensionality

Kx = gpuArray(Kx);
Ky = gpuArray(Ky);

n_modes = (n_max+1)^2;

d2k = (Kx(1,2)-Kx(1,1))*(Ky(2,1)-Ky(1,1));

% first basis function is the aperture function
poly_basis = gpuArray(zeros([size(Kx),n_modes]));
poly_coeff = zeros([n_max+1,n_max+1,n_modes]);

% projection coefficients for GS orthonormalization
proj_coeff = zeros([n_modes,1,n_modes]);  

mode = 0;
for n = 0:n_max
    for m = 0:n_max
        mode = mode+1;
        
        % moments for candidate polynomial
        ikxn = (1i*Kx).^n;
        ikym = (1i*Ky).^m;
        
        % candidate polynomial
        p_nm = ikxn .* ikym;

        % calculate projection coefficients for all previous basis
        % polynomials
        proj_coeff(mode,1,:) = sum(d2k * conj(poly_basis) .* p_nm .* abs(A).^2, [1,2]); 
        
        % remove projections (for all but the 0th mode)
        p_nm = p_nm - sum( proj_coeff(mode,1,:) .* poly_basis,3);
        
        % get polynomial expansion coefficients for basis polynomial
        poly_coeff(:,:,mode) = -1*sum(proj_coeff(mode,1,:) .* poly_coeff, 3);
        poly_coeff(n+1,m+1,mode) = 1;
        
        % normalize polynomial
        N = sum(d2k * abs(p_nm).^2 .* abs(A).^2,'all');
        p_nm = p_nm /sqrt(N); 
        poly_coeff(:,:,mode) = poly_coeff(:,:,mode)/sqrt(N);
                
        % add a new basis polynomial to stack
        poly_basis(:,:,mode) = p_nm;       
    
    end
end

poly_basis = gather(poly_basis);

end

function [GS_basis_pos, GS_basis_mom] = GSPos(A,n_max,ap_dim,ip_dim,rel_ap)

% number of modes
n_modes = (n_max+1)^2;

% GS basis in momentum space representation
GS_basis_mom = poly_basis .* A;


% GS basis in position space representation via FFT
GS_basis_pos = gpuArray(zeros(ip_dim,ip_dim,n_modes));

for mode = 1:n_modes 
    new_dim = round(rel_ap*ap_dim*ip_dim / (2*1.22) );
    % new_dim = round((subap2ap)*ap_dim*ip_dim/(2*1.22)); new_dim = new_dim + (mod(new_dim,2)+ mod(ip_dim,2));
    tilde_phi = padarray(GS_basis_mom(:,:,mode),ceil((new_dim-ap_dim)/2 * [1,1]),0,'both'); 
    phi = fftshift(ifft2(ifftshift(tilde_phi)));
    phi = CenteredRectCrop(phi,ip_dim,ip_dim)*(new_dim/ip_dim)^2;
    GS_basis_pos(:,:,mode) = phi;
end


% gather the gpu arrays
GS_basis_mom = gather(GS_basis_mom);
GS_basis_pos = gather(GS_basis_pos);




end
