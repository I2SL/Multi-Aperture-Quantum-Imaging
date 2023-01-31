function [err,src_coords,est_coords,aper_coords,mode_counts,max_order,rl] = ...
        MultiAperture_ConstellationLocalization(basis,ap_num,src_num,rl_frac,n_pho)
    
    % static measurement estimation of point-source constellations with multi-aperture systems
    

close all
addpath('utils/')


% imaging system parameters
subap_radius = 1.5; % radius of sub-apertures [units of length]

% make the sub-aperture length the reference unit
ref_unit = subap_radius;

% rescale everything to the reference unit
subap_radius = subap_radius/ref_unit;   % radius of sub-apertures  in reference unit: [u]
A_sub = pi*subap_radius^2;              % subaperture collection area [u^2]

% multi-aperture parameters
aper_coords = 4*subap_radius*genNgon(ap_num,0);      % aperture coordinates [u]
A_tot = ap_num * A_sub;                              % total collection area of the multi-aperture system [u^2]
B = pdist(aper_coords);                              % baseline lengths [u]
max_B = max([1,B]);                                  % max baseline [u]

% rayleigh length
rl_1ap = 2*pi*1.22;
rl = rl_1ap/max_B;                                    % rayleigh length [rad/u] 

% image plane discretization  
ip_dim = 101;                                         % image plane pixels along each axis
[X,Y] = meshgrid(rl * linspace(-.5,.5,ip_dim));       % image plane coordinates [rad/u]

% aperture plane discretization
subap_sampling = 301;     % number of k-space samples across each circular hard sub-aperture (total samples is subap_sampling^2)
[Kx,Ky,d2k] = ApertureKxKy(aper_coords,subap_sampling); 

% source distribution
src_coords = rl*rl_frac*genNgon(src_num,1);     % source coordinates [rad/u]
n_sources = size(src_coords,1);                 % number of sources in the scene
s_b = ones(n_sources,1)/n_sources;              % relative source brightnesses
s_x = src_coords(:,1); s_y = src_coords(:,2); 

% mode parameters
max_order = 3;                                      % max number of n-index modes to truncate Hilbert Space at

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
        
        % visualize the modes
        Visualize_GramSchmidt(nj,mj, GS_basis_pos)
        

    case 'Zernike'
        
        % indices
        [nj,mj,vj] = Indices_MixedAperture(max_order,ap_num);
        
        % number of modes
        num_modes = numel(nj);
        
        % Create Mixed-Aperture Zernike Basis
        U = dftmtx(ap_num)/sqrt(ap_num);   % a unitary matrix for mixing aperture-local modes
        
        % probability function handle
        prob_fn = @(xq,yq) ModalProb_MixedAperture([xq,yq],nj,mj,vj,U,aper_coords,A_tot);
        
        % visualize the modes
        Visualize_MixedAperture(nj,mj,vj,X,Y,U,aper_coords);
        

    
    case 'Direct-Detection'
        
        % number of modes
        num_modes = numel(X);
        
        % probability function
        prob_fn = @(xq,yq) ModalProb_DirectImaging([xq,yq],X,Y,aper_coords);
        
        % visualize the PSF
        Visualize_PSF2(3*X,3*Y,rl,aper_coords);        
end

% get modal probabilities for the given source distribution
p_s = sum(s_b .* prob_fn(s_x,s_y),1);

% simulate the measurement
[measurement, mode_counts] = simulateMeasurement(n_pho, p_s);

figure
scatter(Kx,Ky,'filled');            hold on;
scatter(0,0,10,'filled','black');   hold off;
axis 'equal'
title('Aperture Configuration')
xlabel('$k_x \, [\delta]$','interpreter','latex')
ylabel('$k_y \, [\delta]$','interpreter','latex')
axis off


figs(1) = figure;
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
        %xticks([1,ceil(101/2),101])
        %xticklabels({'-1/2','0','rl/2'})
        ylabel('y [rl]')
        %yticks([1,ceil(101/2),101])
        %yticklabels({'-1/2','0','rl/2'})
        axis square
    
end


figs(2) = figure;
scatter(s_x/rl,s_y/rl,'filled','black')
xlim([-.5,.5])
xticks(linspace(-.5,.5,7));
xlabel('x (rl)')
ylim([-.5,.5])
yticks(linspace(-.5,.5,7));
ylabel('y (rl)')
hold on

% find MLE
[s_b_mle, s_x_mle, s_y_mle, count] = EM(measurement,n_sources,prob_fn,X,Y,rl);
est_coords = [s_x_mle,s_y_mle];

% compute the localization error
err = LocalizationError(src_coords, est_coords, rl);
                    
% plot the final constellation estimate
hold on
scatter(s_x_mle/rl,s_y_mle/rl,'filled','red')
hold off
axis square
title({'Expectation Maximization',[num2str(ap_num),'-aperture ', basis]})
names = cell(1,count+3);
names(:) = {''}; names(1) = {'Ground Truth'}; names(end) = {'EM Estimate'};
legend(names)


% save outputs
save_dir = fullfile(basis,[num2str(ap_num),'-aperture'],[num2str(src_num),'-src'],[num2str(rl_frac),'rl']);
mkdir(save_dir)

meas_file = fullfile(save_dir,'measurement.png');
est_file = fullfile(save_dir,'estimate.png');
fig_file = fullfile(save_dir,'figures.fig');
pst_cnt = 1;
while isfile(meas_file) || isfile(est_file) || isfile(fig_file)
    meas_file = fullfile(save_dir,['measurement',num2str(pst_cnt),'.png']);
    est_file =  fullfile(save_dir,['estimate',num2str(pst_cnt),'.png']);
    fig_file = fullfile(save_dir,['figures',num2str(pst_cnt),'.fig']);
    pst_cnt = pst_cnt + 1;
end

% save figures
savefig(figs,fig_file)

% save images
saveas(figs(1),meas_file);
saveas(figs(2),est_file);

end


%% Estimation Functions
function [s_b, s_x, s_y, count] = EM(measurement,num_sources,prob_fn,X,Y,rl)
% runs expectation maximization to determine source coordinates and source
% brightnesses from a measurement.
%
% measurement       : 
% num_sources       : How many sources involved in the scene
% prob_fn           : modal photon detection PMF given a source located at (x,y)
% X
% Y
% rl

% initialize source positions
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
n_em = 20; % max number of EM iterations

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

function [measurement,mode_count] = simulateMeasurement(n_pho)
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
    % of all sources
    
    % calculate the total error distance between all possible pairings of
    % the ground truth sources and the estimated sources
    num_sources = size(xy_src,1);
    P = perms(1:num_sources);   % matrix of all permutations of the source indices
    
    cum_dist = zeros(size(P,1),1);
    
    for j = 1:size(P,1)
         index_permutation = P(j,:);
         xy_delta = xy_src - xy_est(index_permutation,:);
         deltas = vecnorm(xy_delta,2,2); % a vector of Euclidean distances between each source and its estimate
         cum_dist(j) = sum(deltas);
    end
       
    
    % get the minimum distance as the error
    err = min(cum_dist);
    
    % normalize the error by the rayleigh length and the number of sources
    err = err / num_sources / rl; 

end


