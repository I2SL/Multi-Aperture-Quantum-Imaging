function [est_scene,mode_counts,rl,err] = ...
        MultiAperture_ConstellationLocalization(...
        n_pho,...                   % mean photon number                   [integer]
        max_order,...               % max modal order                      [integer]
        basis,...                   % basis                                [string] ['Gram-Schmidt','Zernike', 'Direct-Detection']
        subap_radius,...            % sub-aperture radius                  [double] [units : length]
        aper_coords,...             % aperture position                    [Mx2]    [units : length] --> col_1 = kx, col_2 = ky   
        scene,...                   % input scene                          [Nx3] scene(:,1:2)->source coordinates [in fractions of rayleigh], scene(:,3)->relative source brightnesses
        subap_samp,...              % sub-aperure samples along each axis  --> Bounding box for each aperture has dimensions [subap_samp,subap_samp]
        img_samp,...                % image-plane samples along each axis  --> Image plane for source position estimates has dimensions [img_samp,img_samp]
        EM_max,...                  % max number of EM iterations          [integer]
        visualize...               % visualization trigger                [boolean]
)
% static measurement estimation of point-source constellations with multi-aperture systems
    
% make the sub-aperture length the reference unit
ref_unit = subap_radius;                % [u]

% (non-dimensionalize) rescale aperture-plane coordinates to the reference unit
subap_radius = subap_radius/ref_unit;   % radius of reference sub-apertures [u]
aper_coords = aper_coords/ref_unit;     % sub-aperture coordinates [u]

% multi-aperture parameters
ap_num = size(aper_coords,1);           % number of sub-apertures
A_sub = pi*subap_radius^2;              % subaperture collection area [u^2]
A_tot = ap_num * A_sub;                 % total collection area of the multi-aperture system [u^2]
B = pdist(aper_coords);                 % baseline lengths [u]
max_B = max([1,B]);                     % max baseline [u]

% rayleigh lengths
rl_sub = 2*pi*1.22;                     % sub-aperture rayleigh length [rad/u]
rl = rl_sub/max_B;                      % effective aperture rayleigh length [rad/u] 

% source distribution
% scene = [s_x,s_y,s_b];
src_coords = rl*scene(:,1:2);            % source coordinates [rad/u]
num_sources = size(src_coords,1);        % number of sources in the scene
s_x = src_coords(:,1); s_y = src_coords(:,2); 
s_b = scene(:,3);                        % relative source brightnesses

% image plane discretization
[X,Y] = meshgrid(rl * linspace(-.5,.5,img_samp));       % image plane coordinates [rad/u]

% aperture plane discretization
[Kx,Ky,d2k] = ApertureKxKy(aper_coords,subap_samp);     % Kx [u], Ky [u], d2k [u^2]

% make sure min source separation is greater than the resolution of the image plane
min_sep = min(pdist(src_coords));
dx = X(1,2) - X(1,1);
if min_sep < dx
    warning('Image plane discretization is coarser than minimum source separation')
end


% setup basis
switch basis
    case 'Gram-Schmidt'
        
        % indices
        [nj,mj] = Indices_GramSchmidt(max_order);
        
        % number of modes
        num_modes = numel(nj);
        
        % Create Gram-Schmidt basis
        GS_basis_mom = genGramSchmidtBasis_mom(max_order,Kx,Ky,d2k);                 % basis functions in momentum space
        GS_basis_pos = Basis_GramSchmidt_mom([X(:),Y(:)],Kx,Ky,d2k,GS_basis_mom);    % basis functions in position space
        
        % probability function
        GS_basis_pos = reshape(GS_basis_pos,[size(X),num_modes]);                    % 2D basis matrix stack
        prob_fn = @(xq,yq) ModalProb_GramSchmidt_pos([xq,yq],X,Y,GS_basis_pos,A_tot);
                

    case 'Zernike'
        
        % indices
        [nj,mj,vj] = Indices_MixedAperture(max_order,ap_num);
        
        % number of modes
        num_modes = numel(nj);
        
        % Create Mixed-Aperture Zernike Basis
        U = dftmtx(ap_num)/sqrt(ap_num);   % a unitary matrix for mixing aperture-local modes
        
        % probability function handle
        prob_fn = @(xq,yq) ModalProb_MixedAperture([xq,yq],nj,mj,vj,U,aper_coords,A_tot);

    case 'Direct-Detection'
        
        % number of modes
        num_modes = numel(X);
        
        % probability function
        prob_fn = @(xq,yq) ModalProb_DirectImaging([xq,yq],X,Y,aper_coords);
                
end

% get modal probabilities for the given source distribution
p = sum(s_b .* prob_fn(s_x,s_y),1);

% simulate the measurement
[measurement, mode_counts] = simulateMeasurement(n_pho, p);


% find MLE of scene parameters given the measurement
[s_b_trc, s_x_trc, s_y_trc, count] = EM(measurement,num_sources,prob_fn,X,Y,rl,EM_max);
% intermediate scene parameter estimates
s_b_im = s_b_trc(:,1:count-1); s_x_im = s_x_trc(:,1:count-1); s_y_im = s_y_trc(:,1:count-1);
% final scene parameter estimates
s_b_mle = s_b_trc(:,end); s_x_mle = s_x_trc(:,end); s_y_mle = s_y_trc(:,end);

est_coords = [s_x_mle,s_y_mle];
est_brites = s_b_mle;
est_scene = [est_coords/rl, est_brites];

% compute the localization error
err = LocalizationError(src_coords, est_coords, rl);

% visualize figures
if visualize
    
	% APERTURE
    figs(1) = figure;
    scatter(Kx,Ky,'filled','blue');            hold on;
    scatter(0,0,10,'filled','black');   hold off;
    axis 'equal'
    title('Aperture Configuration')
    xlabel('$k_x \, [\delta]$','interpreter','latex')
    ylabel('$k_y \, [\delta]$','interpreter','latex')
    
    % MODES
    figs(2) = figure;
    switch basis
        case 'Gram-Schmidt'
            % visualize the GS modes
            Visualize_GramSchmidt(nj,mj, GS_basis_pos);
            
        case 'Zernike'
            % visualize the Mixed Zernike modes
            Visualize_MixedAperture(nj,mj,vj,X,Y,U,aper_coords);
        
        case 'Direct-Detection'
            % visualize the PSF
            Visualize_PSF2(3*X,3*Y,rl,aper_coords);
    end
           
    
    % MEASUREMENT
    figs(3) = figure;
    switch basis
        case 'Gram-Schmidt'
            stem(mode_counts);
            title({'Photon Counting Measurement','Gram-Schmidt Basis',['Total Photons: ',num2str(sum(mode_counts))]});
            xlabel('mode index')
            ylabel('# photons')

            n_labels = arrayfun(@num2str,nj,'UniformOutput', 0);
            m_labels = arrayfun(@num2str,mj,'UniformOutput', 0);
            index_labels = strcat(n_labels,repmat({','},[1,num_modes]),m_labels);
            xticks(1:num_modes)
            xticklabels(index_labels)


        case 'Zernike'
            stem(mode_counts);
            title({'Photon Counting Measurement','FT Zernike Basis',['Total Photons: ',num2str(sum(mode_counts))]});
            xlabel('mode index')
            ylabel('# photons')
            n_labels = arrayfun(@num2str,nj,'UniformOutput', 0);
            m_labels = arrayfun(@num2str,mj,'UniformOutput', 0);
            v_labels = arrayfun(@num2str,vj,'UniformOutput', 0);
            index_labels = strcat(n_labels,repmat({','},[1,num_modes]),m_labels,repmat({','},[1,num_modes]),v_labels);
            xticks(1:num_modes)
            xticklabels(index_labels)


        case 'Direct-Detection'
            DD_photons = reshape(mode_counts, ip_dim*[1,1]);
            imagesc([min(X(:)),max(X(:))]/rl,[min(Y(:)),max(Y(:))]/rl,DD_photons)
            colorbar
            title({'Direct Detection Measurement',['Total Photons: ',num2str(sum(mode_counts))]});
            xlabel('x [rl]')
            ylabel('y [rl]')
            axis square

    end  

    % ESTIMATE
    figs(4) = figure;
    % plot ground truth
    scatter(s_x/rl,s_y/rl,'filled','black'); hold on;
    % plot the intermediate estimates
    if count > 1 
        scatter(s_x_im/rl,s_y_im/rl)
    end
    % plot the final constellation estimate
    scatter(s_x_mle/rl,s_y_mle/rl,'filled','red')
    hold off
    title({'Expectation Maximization',[num2str(ap_num),'-aperture ', basis]})
    xlim([min(X(:)),max(X(:))]); xlabel('x (rl)');
    ylim([min(Y(:)),max(Y(:))]); ylabel('y (rl)');
    xticks(linspace(min(X(:)),max(X(:)),7));
    yticks(linspace(min(Y(:)),max(Y(:)),7));
    names = cell(1,count+1);
    names(:) = {''}; names(1) = {'Ground Truth'}; names(end) = {'EM Estimate'};
    legend(names)
    axis square   
    
end

end


%% Estimation Functions
function [s_b_trc, s_x_trc, s_y_trc, count] = EM(measurement,num_sources,prob_fn,X,Y,rl,n_em_max)
% runs expectation maximization to determine source coordinates and source
% brightnesses from a measurement.
%
% measurement       : 
% num_sources       : How many sources involved in the scene
% prob_fn           : modal photon detection PMF given a source located at (x,y)
% X                 :
% Y                 :
% rl                :
% n_em_max          : max number of EM iterations

% initialize source positions

% random sub-rayleigh source position initialization
s_x = rl*(rand(num_sources,1)-.5);
s_y = rl*(rand(num_sources,1)-.5);

% radial sub-rayleigh source position initialization
%s_x = rl/4.*cos(2*pi*(0:num_sources-1)/num_sources)';
%s_y = rl/4.*sin(2*pi*(0:num_sources-1)/num_sources)';

% initialize source weights
s_b = ones(num_sources,1)/num_sources;

% traces of the scene parameter estimates accross EM iterations
s_x_trc = zeros(num_sources,1);
s_y_trc = zeros(num_sources,1);
s_b_trc = zeros(num_sources,1);


count = 0;
while ( ~isempty(find(s_x-s_x_trc(:,end), 1)) || ~isempty(find(s_y-s_y_trc(:,end), 1))) && count<=n_em_max  
    
    % increment count
    count = count + 1;
    
    s_x_trc(:,count) = s_x;
    s_y_trc(:,count) = s_y;
    s_b_trc(:,count) = s_b;
    
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
    
    
    Q_2D = zeros([size(X),num_sources]);
    for i = 1:num_sources
        Q_2D(:,:,i) = reshape(Q(i,:),size(X));
    end 
    
    %{
    % debugging
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
    
    %dx = X(1,2) - X(1,1);
    %dy = Y(2,1) - Y(1,1);
    
    % assign source positions 
    s_x = X(ind_sxy);% + dx*randn(numel(ind_sxy),1); %(plus a bit of noise to avoid degeneracies) 
    s_y = Y(ind_sxy);% + dy*randn(numel(ind_sxy),1); %(plus a bit of noise to avoid degeneracies) 
    
end


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

function [measurement,mode_count] = simulateMeasurement(n_pho,p)
    N = poissrnd(n_pho); % number of photons collected
    
    % randomly assign modal bin to each photon according to PMF
    modes = 1:numel(p);
    measurement = datasample(modes,N,'Weights',p);
    
    % count photons in modal bins
    [gc,gs] = groupcounts(measurement');
    mode_count = zeros(numel(p),1);
    mode_count(gs) = gc;
end

function err = LocalizationError(xy_src,xy_est,rl)
    % Computes the average localization error (in fractions of a rayleigh)
    % accross all source position estimates
    
    % calculate the total error distance between all possible pairings of
    % the ground truth sources and the estimated sources
    num_sources = size(xy_src,1);
    P = perms(1:num_sources);   % matrix of all permutations of the source indices
    
    sum_dist = zeros(size(P,1),1);
    
    for j = 1:size(P,1)
         index_permutation = P(j,:);
         xy_delta = xy_src - xy_est(index_permutation,:);
         deltas = vecnorm(xy_delta,2,2); % a vector of Euclidean distances between each source and its estimate
         sum_dist(j) = sum(deltas);
    end
       
    
    % select the minimum cumulative distance as the true error.
    % this assumes that each estimated source location is matched to the
    % nearest ground truth source
    err = min(sum_dist);
       
    % get the average error per source by dividing the number of sources.
    err = err / num_sources;
    
    % normalize the average error by the rayleigh length
    % to get the average error as a fraction of the rayleigh length
    err = err / rl;
    
end


