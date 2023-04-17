function [est_scene,mode_counts,rl,err] = ...
        MultiAperture_ConstellationLocalization(...
        n_pho,...                   % mean photon number                   [integer]
        max_order,...               % max modal order                      [integer]
        basis,...                   % basis                                [string] ['Gram-Schmidt','Zernike', 'Direct-Detection']
        aperture,...                % aperture coordinates and radii       [Mx3] aperture(:,1:2) --> centroid coodrdinates of sub-apertures, aperture(:,3) radius of each sub-aperture 
        scene,...                   % input scene                          [Nx3] scene(:,1:2)->source coordinates [in fractions of rayleigh], scene(:,3)->relative source brightnesses
        mom_samp,...                % aperture plane sampling density (momentum space)  [integer] [units : samples/length]
        pos_samp,...                % image-plane sampling density (position space)     [integer] [units : samples/rayleigh] Image plane for source position estimates has dimensions [pos_samp,pos_samp]
        EM_max,...                  % max number of EM iterations          [integer]
        brite_flag,...              % estimate brightness trigger          [boolean]
        visualize...                % visualization trigger                [boolean]
)
% static measurement estimation of point-source constellations with multi-aperture systems

% multi-aperture parameters
ap_num = size(aperture,1);
aper_coords = aperture(:,1:2);
aper_rads = aperture(:,3); 

% get the effective aperture diameter
if ap_num>1
    B = squareform(pdist(aper_coords));             % sub-aperture pairwise centroid distances
    D = aper_rads + aper_rads';                     % sub-aperture pairwise radius sums
    assert(all(triu(B,1) >= triu(D,1),'all'));      % check if any apertures overlap.
    
    % set the effective aperture diameter to the minimum enclosing circle diameter.
    cm_coords = aper_coords - mean(aper_coords,1);                          % get centered coordinates (with respect to centroid of centroids -- still need to prove with this works)
    tangent_pts = cm_coords + cm_coords.*aper_rads./vecnorm(cm_coords,2,2); % candidate tangent points where the enclosing circle might touch
    [kx_c,ky_c,R_eff] = SmallestEnclosingCircle(tangent_pts(:,1)',tangent_pts(:,2)');
    D_eff = 2*R_eff;
else
    R_eff = aper_rads(1);
    D_eff = 2*R_eff;                         % set the effective aperture diameter to that of the input aperture
end

% get the rayleigh length of the system
% first zero of Besselj(0,x) should be set so that a source shifted to rl is zero.
bsj1_zero = pi*1.2197;      %= 3.8317 is the first zero of besselj_1
rl = 2*pi * 1.2197/D_eff;   % rayleigh length (in rads/length) - the rayleigh length is the radial distance (in image space) from the origin to the first zero of the besselj(1) function rl = bsj1_zero/R_eff;   

% source distribution
% scene = [s_x,s_y,s_b];
src_coords = rl*scene(:,1:2);            % source coordinates [rad/u]
num_sources = size(src_coords,1);        % number of sources in the scene
s_x = src_coords(:,1); s_y = src_coords(:,2); 
s_b = scene(:,3);                        % relative source brightnesses

% image plane discretization
[X,Y] = meshgrid(rl * linspace(-.5,.5,pos_samp));       % image plane coordinates [rad/u]

% aperture plane discretization
[Kx,Ky,d2k] = ApertureKxKy(aper_coords,aper_rads,mom_samp);     % Kx [u], Ky [u], d2k [u^2]

% make sure min source separation is greater than the resolution of the image plane
min_sep = min(pdist(src_coords));
dx = X(1,2) - X(1,1);
if min_sep < dx
    warning('Image plane discretization is coarser than minimum source separation')
end


% setup basis
switch basis
    case 'Gram-Schmidt'
        
        % assign A_tot to be the area of the discretized aperture
        A_tot = numel(Kx)*d2k;
        
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
        
        % assign A_tot to be the true aperture area
        A_tot = sum(pi*aper_rads.^2);
        
        % indices
        [nj,mj,vj] = Indices_MixedAperture(max_order,ap_num);
        
        % number of modes
        num_modes = numel(nj);
        
        % Create Mixed-Aperture Zernike Basis
        U = dftmtx(ap_num)/sqrt(ap_num);   % a unitary matrix for mixing aperture-local modes
        
        % probability function handle
        %prob_fn = @(xq,yq) ModalProb_MixedAperture([xq,yq],nj,mj,vj,U,aper_coords,A_tot);
        prob_fn = @(xq,yq) ModalProb_MixedAperture([xq,yq],nj,mj,vj,U,aperture);

    case 'Direct-Detection'
        
        % number of modes
        num_modes = numel(X);
        
        % probability function
        prob_fn = @(xq,yq) ModalProb_DirectImaging([xq,yq],X,Y,aperture);
                
end

% get modal probabilities for the given source distribution
p = sum(s_b .* prob_fn(s_x,s_y),1);

% simulate the measurement
isPoiss = 1; % add Poisson arrival statistics to the measurment
[pho_xy_id, mode_counts] = simulateMeasurement(n_pho, p, isPoiss);

% find MLE of scene parameters given the measurement
if strcmp(basis,'Direct-Detection')
    [s_b_trc, s_x_trc, s_y_trc, loglike_trc, count] = EM_DD([X(pho_xy_id)',Y(pho_xy_id)'],num_sources,aperture,rl,EM_max,brite_flag);   
    %[s_b_trc, s_x_trc, s_y_trc, loglike_trc, count] = EM_DD(mode_counts,num_sources,src_coords,prob_fn,X,Y,rl,EM_max,brite_flag);
else
    [s_b_trc, s_x_trc, s_y_trc, loglike_trc, count] = EM(mode_counts,num_sources,src_coords,prob_fn,X,Y,rl,EM_max,brite_flag);
end

% intermediate scene parameter estimates
s_b_im = s_b_trc(:,1:count-1); s_x_im = s_x_trc(:,1:count-1); s_y_im = s_y_trc(:,1:count-1);
% final scene parameter estimates
s_b_mle = s_b_trc(:,end); s_x_mle = s_x_trc(:,end); s_y_mle = s_y_trc(:,end);

% estimated scene
est_coords = [s_x_mle,s_y_mle];
est_brites = s_b_mle;
est_scene = [est_coords/rl, est_brites];

% compute the localization error at each iteration
err_trc = arrayfun(@(k) LocalizationError(src_coords, [s_x_trc(:,k),s_y_trc(:,k)]),1:count);

% final error
err = err_trc(end);

% visualize figures
if visualize
    
	% APERTURE
    figs(1) = figure;
    VisualizeAperture(aperture);
    
    % PSF
    figs(2) = figure;
    X_psf = 4*X; % X range to display psf over
    Y_psf = 4*Y; % Y range to display psf over
    PSF = MultiAperturePSF([X_psf(:),Y_psf(:)],aperture);  % PSF
    PSF2 = reshape(abs(PSF).^2,size(X_psf));               % modulus squared of PSF
    surf(X_psf/rl,Y_psf/rl,PSF2)
    xlabel('x [rl]')
    ylabel('y [rl]')
    title('$|PSF|^2$','interpreter','latex')
    
    % MODES
    switch basis
        case 'Gram-Schmidt'
            % visualize the GS modes
            Visualize_GramSchmidt(nj,mj,X,Y,rl,GS_basis_pos);
            
        case 'Zernike'
            % visualize the Mixed Zernike modes
            %Visualize_MixedAperture(nj,mj,vj,X,Y,rl,U,aper_coords);
            Visualize_MixedAperture(nj,mj,vj,X,Y,rl,U,aperture);
        
        case 'Direct-Detection'
            % visualize the PSF
            Visualize_PSF2(3*X,3*Y,rl,aperture);
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
            DD_photons = reshape(mode_counts, pos_samp*[1,1]);
            imagesc([min(X(:)),max(X(:))]/rl,[min(Y(:)),max(Y(:))]/rl,DD_photons)
            colorbar
            title({'Direct Detection Measurement',['Total Photons: ',num2str(sum(mode_counts))]});
            xlabel('x [rl]')
            ylabel('y [rl]')
            axis square

    end  
    
    
    
    % LOG-LIKELIHOOD VS LOCALIZATION ERROR
    figs(4) = figure;
    subplot(2,1,1)
    plot(loglike_trc)
    title('Likelihood Convergence')
    xlabel('Iteration')
    ylabel('Log Likelihood')
    subplot(2,1,2)
    plot(err_trc/min_sep)
    title('Localization Error Convergence')
    xlabel('Iteration')
    ylabel('Fractional Localization Error')
    
    

    % ESTIMATE
    figs(5) = figure;
    % plot ground truth
    scatter(s_x/rl,s_y/rl,'filled','black'); hold on;
    % plot the intermediate estimates
    if count > 1 
        scatter(s_x_im/rl,s_y_im/rl,5)
    end
    % plot the final constellation estimate
    scatter(s_x_mle/rl,s_y_mle/rl,'red','square')
    hold off
    title({'Expectation Maximization',basis}) 
    xlabel('x (rl)');
    ylabel('y (rl)');
    xlim([min(X(:))/rl,max(X(:))/rl]); 
    ylim([min(Y(:))/rl,max(Y(:))/rl]); 
    xticks(linspace(min(X(:))/rl,max(X(:))/rl,5));
    yticks(linspace(min(Y(:))/rl,max(Y(:))/rl,5));
    names = cell(1,count+1);
    names(:) = {''}; names(1) = {'Ground Truth'}; names(end) = {'EM Estimate'};
    legend(names)
    axis square
    
end

end