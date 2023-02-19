%% A BUNCH OF MISCELLANEOUS FUNCTIONS THAT WERE USEFUL AT ONE POINT AND NOW CAN PROBABLY BE DELETED

x=0;

%{
% aperture plane discretization
ap_dim = 101;
[aperture,Kx,Ky] = ApertureConfig(a_kx,a_ky,ap_dim);
rel_ap = subap_radius / Kx(1,end);                      % ratio of sub-aperture radius to aperture plane half-width
%}

 %{
[poly_coeff,GS_basis_mom,GS_basis_pos] = genGramSchmidtBasis(Kx,Ky,aperture,n_max,ap_dim,ip_dim,rel_ap);
% ensure the position space modes are properly normalized
% (an on-axis source should produce probability 1 in the 00 mode)
GS_normalization = sqrt(A_tot)*abs(Basis_GramSchmidt(0,0,X,Y,GS_basis_pos(:,:,1)));
GS_basis_pos = GS_basis_pos / GS_normalization;
% visualize the modes
[nj,mj] = Indices_GramSchmidt(n_max);
VisualizeModes_GramSchmidt(nj,mj, GS_basis_pos)
% modal probability function
prob_fn = @(x,y) ModalProb_GramSchmidt(x,y,X,Y,GS_basis_pos,A_tot);
%}


function xy_out = removeClose(xy_in,dist)
    % removes coordinates in the list xy_in that are less than dist from
    % each other

    if size(xy_in,1) == 1
        xy_out = xy_in;
        return
    end
    
    D = pdist(xy_in);
    Z = squareform(D);
    
    
    k_list = 1;

    k_val = find(Z(:,k_list(end)) > dist);
    k_val = setxor(k_val,intersect(k_val,k_list));
    
    k_inv = find(Z(:,k_list(end)) <= dist);  
    
    while ~isempty(k_val)
        
        k_list = [k_list, k_val(1)];
        
        k_val = find(Z(:,k_list(end)) > dist);
        k_val = setxor(k_val,intersect(k_val,k_list));
        k_val = setxor(k_val,intersect(k_val,k_inv));
        
        k_inv = unique([k_inv,find(Z(:,k_list(end)) <= dist)]);
    end
    
    xy_out = xy_in(k_list,:);
end

function xy_out = combineClose(xy_in,dist)
    
    if size(xy_in,1) == 1
        xy_out = xy_in;
        return
    end

    % get candidate clusters
    [idx] = dbscan(xy_in,dist,1);
    clusters = unique(idx);
    
    % find mean of clusters and use as combined candidate point
    num_clusters = numel(clusters);
    xy_out = zeros(num_clusters,2);
    
    for c = 1:num_clusters
        xy_out(c,:) = mean(xy_in(idx == clusters(c),:),1);
    end
end

function out_path = treepath(H,in_path)
% searches all candidate source positions and returns the index list that
% corresponds to a unique choice of source locations. If no such list is 
% possible then the function returns -1.

% H  is a cell array containing the list of candidate source locations for
% each source.

h = H{1};                         % max Q source position indices for current source
if isempty(in_path)
    c = h;
else
    common = intersect(in_path,h);
    c = setxor(h,common);   % valid source position indices
end

if isempty(c)
    out_path = -1;
else
    if length(H) == 1
        out_path = [in_path, c(1)];
    else
        out_path = -1;
        k = 1;
        while sum(out_path == -1) && k <= length(c)
            try_path = [in_path, c(k)];
            out_path = treepath(H(2:end),try_path);
            k = k + 1;
        end
    end
end

end

function U = nDimRotationUnitary(n)
    % n-dimensional rotation matrix in the hyperplane spanned by n1,n2
    % rotation angle alpha
    n1 = zeros(n,1); n1(1) = 1; 
    n2 = zeros(n,1); n2(2) = 1;
    
    angle = pi/4;
    U = eye(n)+(n2*n1' - n1*n2') * sin(angle) + (n1*n1' + n2*n2') * (cos(angle)-1);
    

end


function idx_sxy = getMLESourceIndices(Q_2D)
    
    dim_2D = size(Q_2D(:,:,1));
    
    H = cell(1,size(Q_2D,3));
    
    for i = 1:size(Q_2D,3)
        Q_i = Q_2D(:,:,i); 
        [s_iy,s_ix] =  find( Q_i == max(Q_i(:)));
        hi = sub2ind(dim_2D,s_iy,s_ix);
        H{i} = hi;
    end
    
    
    % get all possible MLE constellation candidates
    D = H;
    [D{:}] = ndgrid(H{:});
    Z = cell2mat(cellfun(@(m)m(:),D,'uni',0));
        
    % choose the candidate constellation where the sources are most spread out
    spreads = zeros(size(Z,1),1);
    for i = 1:numel(spreads)
        [s_y,s_x] = ind2sub(dim_2D,Z(i,:));
        spreads(i) = sum(pdist([s_x',s_y']));
    end
    [~,k] = max(spreads);
    idx_sxy = Z(k,:)';
    
    %{
    % debugging
    figure
    for i = 1:size(Q_2D,3)
        [s_iy,s_ix] = ind2sub(dim_2D,H{i});
        scatter(s_ix/101 - 0.5,s_iy/101 - 0.5) 
        hold on
    end
    
    hold off
    title('Candidate Source Locations')
    xlim([-.5,.5])
    ylim([-.5,.5])
    leg = legend( sprintfc('%g', 1:size(Q_2D,3)) );
    title(leg, 'Source')
    
    
    
    [s_y,s_x] = ind2sub(dim_2D,idx_sxy);
    figure
    scatter(s_x/101 - 0.5,s_y/101 - 0.5)
    xlim([-.5,.5])
    ylim([-.5,.5])
    title('Selected Source Locations')
    
    
    % alternatively we find the candidate constellation using a tree search
    % that chooses the first candidate exhibiting no overlap
    % idx_sxy = treepath(H,[])';
    %}
    
end


function I_out =  CenteredRectCrop(I,rect_r,rect_c)
    % crops the matrix I to a rectangle of dimensions [rect_r, rect_c]
    % that is centered on the matrix.
    
    row_max = size(I,1);
    col_max = size(I,2); 
    
    center = [ceil(row_max/2),ceil(col_max/2)];
    delta = [floor(rect_r/2),floor(rect_c/2)];
    
    I_out = I( (center(1)-delta(1)) : (center(1)+delta(1)) ,...
               (center(2)-delta(2)) : (center(2)+delta(2)));
end



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
