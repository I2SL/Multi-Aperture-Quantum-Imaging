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



%% Miscellaneous functions (to delete)
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