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

figure
scatter(s_x/rl,s_y/rl,'filled','black'); 

% simulate the measurement
[~, mode_counts] = simulateMeasurement(n_pho, p);


% find MLE of scene parameters given the measurement
[s_b_trc, s_x_trc, s_y_trc, count] = EM(mode_counts,num_sources,prob_fn,X,Y,rl,EM_max);
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
    title('Aperture Configuration','interpreter','latex')
    xlabel('$k_x \, [\delta]$','interpreter','latex')
    ylabel('$k_y \, [\delta]$','interpreter','latex')
    
    % MODES
    figs(2) = figure;
    switch basis
        case 'Gram-Schmidt'
            % visualize the GS modes
            Visualize_GramSchmidt(nj,mj,GS_basis_pos);
            
        case 'Zernike'
            % visualize the Mixed Zernike modes
            Visualize_MixedAperture(nj,mj,vj,X,Y,rl,U,aper_coords);
        
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
            DD_photons = reshape(mode_counts, img_samp*[1,1]);
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