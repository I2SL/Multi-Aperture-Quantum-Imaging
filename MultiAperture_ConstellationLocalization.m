function MultiAperture_ConstellationLocalization(method,ap_num,src_num,rl_frac,n_pho)

close all

% multi-aperture static measurement estimation
nu = 1;             % flux attenuation 

% imaging system parameters
subap_radius = 0.5; % radius of sub-apertures [units of length]

% make the sub-aperture length the reference unit
ref_unit = subap_radius;

% rescale everything to the reference unit
subap_radius = subap_radius/ref_unit;  % radius of sub-apertures  [reference unit]

% rayleigh length
rl = 2*pi*1.22/(subap_radius);        % rayleigh length
A_sub = pi*subap_radius^2;              % subaperture collection area

%[a_kx,a_ky] =  genGolay9(D);           % aperture coordinates [reference units]
aper_coords = 4*genNgon(ap_num,0); 
%aper_coords = [0,0];
a_kx = aper_coords(:,1); a_ky = aper_coords(:,2);
n_apertures = size(a_kx,1);             % number of sub-apertures
A_tot = n_apertures * A_sub;            % total collection area of the multi-aperture system

% image plane discretization  
ip_dim = 101;
[X,Y] = meshgrid(rl * linspace(-.5,.5,ip_dim));

% aperture plane discretization
ap_dim = 101;
[aperture,Kx,Ky] = ApertureConfig(a_kx,a_ky,ap_dim);
rel_ap = subap_radius / Kx(1,end);     % ratio of sub-aperture diameter to aperture plane width

% scene parameters 

% source coordinates
src_coords = rl*rl_frac*genNgon(src_num,0);
n_sources = size(src_coords,1);           % number of sources in the scene
s_x = src_coords(:,1); s_y = src_coords(:,2);
s_b = ones(n_sources,1)/n_sources;             % source brightnesses



% mode parameters
n_max = 3;         % max number of radial modes to truncate Hilbert Space with


switch method
    case 'Gram-Schmidt'

        [poly_coff,GS_basis_mom] = GSbasis_2DCartesian(Kx,Ky,aperture,n_max);
        
        for mode = 1:size(GS_basis_mom,3) 
            new_dim = round(rel_ap*ap_dim*ip_dim / (2*1.22) );
            % new_dim = round((subap2ap)*ap_dim*ip_dim/(2*1.22)); new_dim = new_dim + (mod(new_dim,2)+ mod(ip_dim,2));
            tilde_phi = padarray(GS_basis_mom(:,:,mode).*aperture,ceil((new_dim-ap_dim)/2 * [1,1]),0,'both'); 
            phi = fftshift(ifft2(ifftshift(tilde_phi)));
            phi = cropSquare(phi,ip_dim,ip_dim)*(new_dim/ip_dim)^2;
            GS_basis_pos(:,:,mode) = phi;
        end
 
        % rescale the probabilitiy to be commensurate with PSF (shouldn't
        % have to do this but its a temporary hack)
        PSF = reshape(directImagingModalProb(0,0,X,Y,a_kx,a_ky),[ip_dim,ip_dim]);
        scaling = PSF(ceil(ip_dim/2),ceil(ip_dim/2)) / (abs(GS_basis_pos(ceil(ip_dim/2),ceil(ip_dim/2),1))^2);
        GS_basis_pos = sqrt(scaling) * GS_basis_pos;
        
        % visualize the modes
        %[nj,mj] = gramschmidtIndices(n_max);
        %VisualizeModes(nj,mj, GS_basis_pos)
        
        % modal probability function
        prob_fn = @(x,y) gramschmidtModalProb(x,y,X,Y,GS_basis_pos,A_tot);

    case 'Zernike'
        
        [nj,mj,vj] = zernikeIndices(n_max,n_apertures);
        prob_fn = @(x,y) ModalProb(x,y,nj,mj,vj,a_kx,a_ky,A_tot);
    
    case 'Direct-Detection'
    
        prob_fn = @(x,y) directImagingModalProb(x,y,X,Y,a_kx,a_ky);

        
end

% get modal probabilities for the given source distribution
p_s = sum(s_b .* prob_fn(s_x,s_y),1);

% simulate the measurement
[measurement, mode_counts] = simulateMeasurement(n_pho, nu, p_s);

figure
imagesc(flipud(aperture))
colormap('gray')
title('Aperture Configuration')
xlabel('$k_x \, (\delta)$','interpreter','latex')
ylabel('$k_y \, (\delta)$','interpreter','latex')
axis off
axis square


meas_fig = figure(1);
switch method
    case 'Gram-Schmidt'
        stem(mode_counts);
        title({'Photon Counting Measurement','Gram-Schmidt Basis',['Total Photons: ',num2str(sum(mode_counts))]});
        xlabel('mode index')
        ylabel('# photons')
        
    case 'Zernike'
        stem(mode_counts);
        title({'Photon Counting Measurement','FT Zernike Basis',['Total Photons: ',num2str(sum(mode_counts))]});
        xlabel('mode index')
        ylabel('# photons')
    case 'Direct-Detection'
        DD_photons = reshape(mode_counts, ip_dim*[1,1]);
        imagesc(flipud(DD_photons))
        colorbar
        title({'Direct Detection Measurement',['Total Photons: ',num2str(sum(mode_counts))]});
        xlabel('x (rl)')
        xticks([1,ceil(101/2),101])
        xticklabels({'-1/2','0','rl/2'})
        ylabel('y (rl)')
        yticks([1,ceil(101/2),101])
        yticklabels({'-1/2','0','rl/2'})
        axis square
    
end



est_fig = figure(2);
scatter(s_x/rl,s_y/rl,'filled','black')
xlim([-1/2,1/2])
xticks([-0.5,0,0.5]);
xticklabels({'-1/2','0','1/2'})
xlabel('x (rl)')
ylim([-1/2,1/2])
yticks([-0.5,0,0.5]);
yticklabels({'-1/2','0','1/2'})
ylabel('y (rl)')
hold on

% find MLE
[s_b_mle, s_x_mle, s_y_mle, count] = EM(measurement,n_sources,prob_fn,X,Y);

hold on
scatter(s_x_mle/rl,s_y_mle/rl,'filled','red')
hold off
axis square
title({'Expectation Maximization',[num2str(ap_num),'-aperture ', method]})
names = cell(1,count+3);
names(:) = {''}; names(1) = {'Ground Truth'}; names(end) = {'EM Estimate'};
legend(names)


% save outputs
save_dir = fullfile(method,[num2str(ap_num),'-aperture'],[num2str(src_num),'-src'],[num2str(rl_frac),'rl']);
mkdir(save_dir)

meas_file = fullfile(save_dir,'measurement.png');
est_file = fullfile(save_dir,'estimate.png');
pst_cnt = 1;
while isfile(meas_file) || isfile(est_file)
    meas_file = fullfile(save_dir,['measurement',num2str(pst_cnt),'.png']);
    est_file =  fullfile(save_dir,['estimate',num2str(pst_cnt),'.png']);
    pst_cnt = pst_cnt + 1;
end

saveas(meas_fig,meas_file);
saveas(est_fig,est_file);

end

% simulate a measurement
function [measurement,mode_count] = simulateMeasurement(n_pho,nu,p)
    N = poissrnd(n_pho*nu); % number of photons collected
    
    % randomly assign modal bin to each photon according to PMF
    modes = 1:numel(p);
    measurement = datasample(modes,N,'Weights',p);
    
    % count photons in modal bins
    [gc,gs] = groupcounts(measurement');
    mode_count = zeros(numel(p),1);
    mode_count(gs) = gc;
end


function [nj,mj,vj] = zernikeIndices(n_max,n_apertures)
    nj = [];
    mj = [];

    for n = 0:n_max
        for m = -n:2:n
            nj = [nj, n];
            mj = [mj, m];
        end
    end

    n_modes = numel(nj);

    % radial, azimuthal, and aperture index lists
    [vj,jj] = ndgrid(1:n_apertures,1:n_modes);
    nj = nj(jj(:));
    mj = mj(jj(:));
    vj = vj(:)';
    
end

function [nj,mj] = gramschmidtIndices(n_max)
    nj = [];
    mj = [];

    for n = 0:n_max
        for m = 0:n_max
            nj = [nj, n];
            mj = [mj, m];
        end
    end
    %{
    n_modes = numel(nj);

    % radial, azimuthal, and aperture index lists
    [vj,jj] = ndgrid(1:n_apertures,1:n_modes);
    nj = nj(jj(:));
    mj = mj(jj(:));
    vj = vj(:)';
    %}
end

function [A,Kx,Ky] = ApertureConfig(a_kx,a_ky,ap_dim)
    % gives the multi-aperture function
    % a_kx and a_ky are the coordinates of the sub-apertures in fractional
    % amounts of the aperture radius (dimensionless).
    
    kx_delta = max(abs(a_kx))+1;
    ky_delta = max(abs(a_ky))+1;
    delta = max(kx_delta,ky_delta);
    
    % aperture plane coordinate grid
    [Kx,Ky] = meshgrid(linspace(-delta,delta,ap_dim));
    dkx = (Kx(1,2)-Kx(1,1)); dky =(Ky(2,1)-Ky(1,1));
    
    circ = @(u,v) (u.^2 + v.^2 < 1);
    
    % construct aperture
    A = zeros(size(Kx));
    for k = 1:numel(a_kx)
        A = A + circ(Kx-a_kx(k),Ky - a_ky(k));
    end
    
    A = A / sqrt(sum(dkx*dky * abs(A).^2, 'all'));   
end

function p = gramschmidtModalProb(x,y,X,Y,GS_modes,A_tot)
    
    dx = X(1,2) - X(1,1);
    dy = Y(2,1) - Y(1,1);

    n_modes = size(GS_modes,3);
        
    p = zeros([numel(x),n_modes]);
    for i = 1:n_modes
        GS_i = GS_modes(:,:,i);
        p(:,i) = interp2(X,Y, (2*pi / sqrt(A_tot) * abs(GS_i)).^2,x,y) * dx * dy;
    end
    
    p = p ./ sum(p,2);              % normalize
end


function phi_nm = ModalBasis(x,y,n,m,v,a_kx,a_ky)
    % image plane coordinates
    % x - [Nx1] 
    % y - [Nx1]
    % mode indices
    % n - [1xM]  radial order index
    % m - [1xM]  azimuthal order index
    % a - [1xM] aperture index
    % sub-aperture positions
    % a_kx - [Px1]
    % a_ky - [Px1]
    %--------------
    % phi_nm - [NxM]

    if numel(a_kx) == 1
        phase = exp(1i*(a_kx .*x + a_ky .*y));
    else
        phase = exp(1i*(a_kx(v)'.*x + a_ky(v)'.*y));
    end
    [th,r] = cart2pol(x,y);
    phi_nm = FTZernikeBasis(r,th,n,m) .* phase;    

end


function p = directImagingModalProb(s_x,s_y,X,Y,a_kx,a_ky)
    % Defines the probability distribution over the image plane for
    % given source positions s_x, s_y and aperture positions a_kx, a_ky.
    % In direct detection the measurement modes are delta functions 
    % on the image plane.
    
    dx = X(1,2) - X(1,1);
    dy = Y(2,1) - Y(1,1);
    p = zeros(numel(s_x),numel(X));
    
    n_apertures = numel(a_kx);
    nj = zeros(1,n_apertures);
    mj = zeros(1,n_apertures);
    vj = (1:n_apertures)';
    
    for k = 1:numel(s_x)
        p(k,:) = abs(sum(ModalBasis(X(:)-s_x(k),Y(:)-s_y(k),nj,mj,vj,a_kx,a_ky),2)).^2 *dx*dy;    
    end
    
end


function p = ModalProb(x,y,nj,mj,vj,a_kx,a_ky,A_tot)
    % returns a discrete PMF over the 2D modes which constitutes the photon modal detection
    % probabilities for photons emitted by a point source at the position x,y

    % correlation function
    correlation_fn = @(x,y)  2*pi / sqrt(A_tot) * conj(ModalBasis(x,y,nj,mj,vj,a_kx,a_ky));

    p = abs(correlation_fn(x,y)).^2;
    p = [p, max(0,1-sum(p,2))];     % add bucket mode due to truncation
    p = p ./ sum(p,2);              % normalize
end

function VisualizeModes(nj,mj,phi_nm)
    
    nn = unique(nj);
    mm = unique(mj);

    
    figure
    for j = 1:numel(nj)
        
        subplot(numel(nn),numel(mm),j)
        
        imagesc(abs(phi_nm(:,:,j)).^2)
        title(['(',num2str(nj(j)),',',num2str(mj(j)),')'])
        axis square

        xticks([1,ceil(size(phi_nm,2)/2),size(phi_nm,2)]);
        xticklabels([-1/2,0,1/2])
        xlabel('rl')

        yticks([1,ceil(size(phi_nm,1)/2),size(phi_nm,1)]);
        yticklabels([-1/2,0,1/2])
        ylabel('rl')
    end    
end

function z = FTZernikeBasis(r,theta,n,m)
    % returns the evaluation of the function and its linear indices
    %
    % INPUTS
    % n - radial index          (row vector 1xK)
    % m - azimuthal index       (row vector 1xK)
    % r - radial argument       (col vector Dx1)
    % theta - angular argument  (col vector Dx1)
    % (r,theta) are coordinate pairs, (n,m) are index pairs
    %
    % OUTPUTS
    % z - the function output   (matrix DxK)
    % j - the linear index      (row vector 1xK)
    
    %z = radius*FTzRadial(radius*r,n) .* FTzAngle(theta,m);    
    z = FTzRadial(r,n) .* FTzAngle(theta,m);    


end

function u = FTzRadial(r,n)
    
    % sinc-bessel in polar
    J = besselj(repmat(n+1,[size(r,1),1]), repmat(r,[1,size(n,2)]))./ repmat(r,[1,size(n,2)]);
    
    % fill in singularities
    J(r==0,n+1 == 1) = 0.5;
    J(r==0,n+1 > 1) = 0;
    
    % radial function
    u = (-1).^(n./2) .* sqrt((n+1)/pi) .*  J;
end

function v = FTzAngle(theta,m)
    v = zeros(numel(theta),numel(m));
    
    % angular function
    c = cos(abs(m).*theta);
    s = sin(abs(m).*theta);
    
    v(:,m>0) = sqrt(2) * c(:,m>0);
    v(:,m==0)= 1;
    v(:,m<0) = sqrt(2) * s(:,m<0);
end

function coords_xy = genRandSource(rl,n_sources)
    x = rl*(rand(n_sources,1)-0.5);
    y = rl*(rand(n_sources,1)-0.5);
   
    coords_xy = [x,y];
end

function coords_xy = genNgon(n,rand_rotation)
    % generates the vertex coordinates of polygon with n vertices equally spaced
    % around a disk. The vertex separation is defined to be 1.
    if n == 1
        coords_xy = [0,0];
        return
    end
    
    
    r = ones([n,1])/2;
    th = 2*pi*(0:n-1)'/n + pi/2;
    del_th = th(2)-th(1);
    
    if rand_rotation
        th = th + 2*pi*rand(1);
    end
    
    % source positions on a unit disk (radius = 1)
    [x,y] = pol2cart(th,r); 
    coords_xy = [x,y];
    
    % make vertex separation 1 by rescaling the radius    
    if n > 2
        coords_xy = coords_xy / sin(del_th/2);
    end
end

function aper_coords  = genGolay9(D)
    
    R1 = D * 1/3 * ones(3,1); 
    R2 = D * 2/3 * ones(3,1);
    R3 = D * 3/3 * ones(3,1);
    rho = [R1; R2; R3];   % golay-9 radii 

    tri = linspace(0,(2*pi)*2/3,3)';
    a1 = (2*pi/3)*0/3 + tri;
    a2 = (2*pi/3)*1/3 + tri;
    a3 = (2*pi/3)*2/3 + tri;
    phi = [a1; a2; a3];   % golay-9 angles 
    
    [a_kx, a_ky] = pol2cart(phi,rho);
    aper_coords = [a_kx,a_ky];
end



function [s_b, s_x, s_y, count] = EM(measurement,num_sources,prob_fn,X,Y)
% runs expectation maximization to determine source coordinates and source
% brightnesses from a measurement.
%
% measurement       :
% num_sources       : How many sources involved in the scene
% prob_fn           : modal photon detection PMF given a source located at (x,y)

% initialize source positions
rl = X(1,end)-X(1,1);
dx = X(1,2) - X(1,1);
dy = Y(2,1) - Y(1,1);
s_x = rl/3*(rand(num_sources,1)-.5);
s_y = rl/3*(rand(num_sources,1)-.5);

%s_x = rl/4.*cos(2*pi*(0:num_sources-1)/num_sources)';
%s_y = rl/4.*sin(2*pi*(0:num_sources-1)/num_sources)';

% initialize source weights
s_b = ones(num_sources,1)/num_sources;

% holders
s_x_old = zeros(size(s_x));
s_y_old = zeros(size(s_y));
s_b_old = zeros(size(s_b));


count = 0;
n_em = 10; % max number of EM iterations

scatter(s_x/rl,s_y/rl);

while ( ~isempty(find(s_x-s_x_old, 1)) || ~isempty(find(s_y-s_y_old, 1))) && count<n_em  
    
    s_x_old = s_x;
    s_y_old = s_y;
    s_b_old = s_b;
    
    % EXPECTATION STEP
    
    % get modal PMFs for each of the source position estimates 
    p_j_est = prob_fn(s_x,s_y); 
      
    % measurement weights
    p_c = s_b .* p_j_est;
    T = p_c(:,measurement);
    T = T ./ sum(T,1);
    
    % get Q
    P_xy = prob_fn(X(:),Y(:));      % PMF for candidate source locations X,Y
    lnP_jxy = log(P_xy(:,measurement));         % probability of detected photon mode at all candidate source locations
    Q = T * lnP_jxy';
    
    %Q = Q + 1e-2*randn(size(Q)); % add gaussian noise to MLE break ties
    
    Q_2D = zeros([size(X),num_sources]);
    for i = 1:num_sources
        Q_2D(:,:,i) = reshape(Q(i,:),size(X));
    end 
    
    %{
    % visualization of objective function for MLE estimators
    peaked_Q = zeros([size(X),num_sources]);
    for i = 1:num_sources
        ki =  find(Q(i,:) == max(Q(i,:))); % max indices
        peaked_Qi = reshape(Q(i,:),size(X));
        peaked_Qi(ki) = peaked_Qi(ki) + 1e4;
        peaked_Q(:,:,i) = peaked_Qi; 
    end 
    %}
    
    % MAXIMIZATION STEP
    %s_b = sum(T,2)/size(T,2); % update source weights
    
    % get max likelihood location indices
    % the following ensures distinct source locations are chosen
    ind_sxy = getMLESourceIndices(Q_2D);
    
    % assign source positions (plus a bit of noise to avoid degeneracies) 
    s_x = X(ind_sxy);% + dx*randn(numel(ind_sxy),1);
    s_y = Y(ind_sxy);% + dy*randn(numel(ind_sxy),1);

    % increment count
    count = count + 1;
    
    % plot intermediate estimate
    scatter(s_x/rl,s_y/rl);
end

hold off
end

function idx_sxy = getMLESourceIndices(Q_2D)
    
    dim_2D = size(Q_2D(:,:,1));
    
    for i = 1:size(Q_2D,3)
        Q_i = Q_2D(:,:,i); 
        [s_iy,s_ix] =  find( Q_i == max(Q_i(:)));
        %s_ixy = round(removeClose([s_ix,s_iy],5)); 
        %s_ix = s_ixy(:,1); s_iy = s_ixy(:,2);
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

    
    % alternatively we find the candidate constellation using a tree search
    % that chooses the first candidate exhibiting no overlap
    % idx_sxy = treepath(H,[])';
    
    %{
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
    %}
    
end

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

function I_out =  cropSquare(I,sq_r,sq_c)
    
    row_max = size(I,1);
    col_max = size(I,2); 
    
    center = [ceil(row_max/2),ceil(col_max/2)];
    delta = [floor(sq_r/2),floor(sq_c/2)];
    
    I_out = I( (center(1)-delta(1)) : (center(1)+delta(1)) ,...
               (center(2)-delta(2)) : (center(2)+delta(2)));
end