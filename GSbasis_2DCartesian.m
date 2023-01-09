%% Generates an arbitrary 2D PAD basis for any aperture function
function [poly_coeff,poly_basis] = GSbasis_2DCartesian(kx,ky,A,n_max)
% AF - aperture function. A function handle describing the aperture
% n_max - max polynomial order

kx = gpuArray(kx);
ky = gpuArray(ky);

n_modes = (n_max+1)^2;

d2k = (kx(1,2)-kx(1,1))*(ky(2,1)-ky(1,1));

% first basis function is the aperture function
poly_basis = gpuArray(zeros([size(kx),n_modes]));
poly_coeff = zeros([n_max+1,n_max+1,n_modes]);

% projection coefficients for GS orthonormalization
proj_coeff = zeros([n_modes,1,n_modes]);  

mode = 0;
for n = 0:n_max
    for m = 0:n_max
        mode = mode+1;
        
        % moments for candidate polynomial
        ikxn = (1i*kx).^n;
        ikym = (1i*ky).^m;
        
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

