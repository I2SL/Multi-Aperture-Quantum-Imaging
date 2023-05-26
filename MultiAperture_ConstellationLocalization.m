function [est_scene,mode_counts,rl,err] = ...
        MultiAperture_ConstellationLocalization(...
        scene,...                   % input scene                          [Nx3] scene(:,1:2)->source coordinates [in fractions of rayleigh], scene(:,3)->relative source brightnesses    
        aperture,...                % aperture coordinates and radii       [Mx3] aperture(:,1:2) --> centroid coodrdinates of sub-apertures, aperture(:,3) radius of each sub-aperture 
        n_pho,...                   % mean photon number                   [integer]
        varargin)


%%%%%%%%%% PARSER %%%%%%%%%%%

% defaults
default_basis = 'Gram-Schmidt';
default_max_order = 5;
default_mom_samp = 167;
default_pos_samp = 129;
default_EM_max = 30;
default_dark_lambda = 0;
default_phase_sigma = 0;
default_brite_flag = 0;
default_exoplanet_flag = 0;
default_visualize = 1;

% initialize parser
P = inputParser;

addRequired(P,'scene')
addRequired(P,'aperture')
addRequired(P,'n_pho')

addOptional(P,'max_order',default_max_order)
addOptional(P,'basis',default_basis)
addOptional(P,'mom_samp', default_mom_samp)
addOptional(P,'pos_samp', default_pos_samp)
addOptional(P,'EM_max', default_EM_max)
addOptional(P,'dark_lambda',default_dark_lambda)
addOptional(P,'phase_sigma',default_phase_sigma)
addOptional(P,'brite_flag', default_brite_flag)
addOptional(P,'exoplanet_flag', default_exoplanet_flag)
addOptional(P,'visualize', default_visualize)

parse(P,scene,aperture,n_pho,varargin{:});


scene =         P.Results.scene;
aperture =      P.Results.aperture;
n_pho =         P.Results.n_pho;
max_order =     P.Results.max_order;
basis =         P.Results.basis;
mom_samp =      P.Results.mom_samp;
pos_samp =      P.Results.pos_samp;
EM_max =        P.Results.EM_max;
dark_lambda =   P.Results.dark_lambda;
phase_sigma =   P.Results.phase_sigma;
brite_flag =    P.Results.brite_flag;
exoplanet_flag = P.Results.exoplanet_flag;
visualize =     P.Results.visualize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% static measurement estimation of point-source constellations with multi-aperture systems


% Checks
if phase_sigma > 0
    assert(strcmp(basis,'Zernike'))
end

if exoplanet_flag
    assert(all(scene(1,1:2) == [0,0]))
end


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
src_brites  = scene(:,3);                % relative source brightnesses
num_sources = size(src_coords,1);        % number of sources in the scene
s_x = src_coords(:,1); s_y = src_coords(:,2); 
s_b = scene(:,3);                        % relative source brightnesses

% image plane discretization
[X,Y] = meshgrid(rl * linspace(-.5,.5,pos_samp));       % image plane coordinates [rad/u]

% aperture plane discretization
[Kx,Ky,d2k] = ApertureKxKy(aper_coords,aper_rads,mom_samp);     % Kx [u], Ky [u], d2k [u^2]

% determine the minimum separation between the sources (if one source, set
% to 1)
if size(scene,1) == 1
    min_sep = 1;
else
    min_sep = min(pdist(src_coords));    
end

% make sure min source separation is greater than the resolution of the image plane
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
        
        prob_fn_measurement = prob_fn;
                

    case 'Zernike'
        
        % assign A_tot to be the true aperture area
        A_tot = sum(pi*aper_rads.^2);
        
        % indices
        [nj,mj,vj] = Indices_MixedAperture(max_order,ap_num);
        
        % number of modes
        num_modes = numel(nj);
        
        % Create Mixed-Aperture Zernike Basis
        U = dftmtx(ap_num)/sqrt(ap_num);   % a unitary matrix for mixing aperture-local modes
        
        % Add fiber delays to the Mixed-Aperture mixing matrix
        U_noisy = U.*exp(1i * 2*pi * normrnd(0,phase_sigma,size(U)));
        
        % probability function handle
        %prob_fn = @(xq,yq) ModalProb_MixedAperture([xq,yq],nj,mj,vj,U,aper_coords,A_tot);
        prob_fn = @(xq,yq) ModalProb_MixedAperture([xq,yq],nj,mj,vj,U,aperture);
        
        prob_fn_measurement = @(xq,yq) ModalProb_MixedAperture([xq,yq],nj,mj,vj,U_noisy,aperture);

    case 'Direct-Detection'
        
        % number of modes
        num_modes = numel(X);
        
        % probability function
        prob_fn = @(xq,yq) ModalProb_DirectImaging([xq,yq],X,Y,aperture);
        
        prob_fn_measurement = prob_fn;
                
end

% get modal probabilities for the given source distribution
p = sum(s_b .* prob_fn_measurement(s_x,s_y),1);



if strcmp(basis,'Direct-Detection')

    % simulate the measurement
    [mode_counts,pho_xy_id] = simulateMeasurement(n_pho, p,...
                                                'isPoiss',1,...
                                                'dark_lambda',dark_lambda );
    
    % find MLE of scene parameters given the measurement
    [s_b_trc, s_x_trc, s_y_trc, loglike_trc, count] = EM_DD([X(pho_xy_id)',Y(pho_xy_id)'],num_sources,aperture,rl,EM_max,brite_flag);   
    %[s_b_trc, s_x_trc, s_y_trc, loglike_trc, count] = EM_DD(mode_counts,num_sources,src_coords,prob_fn,X,Y,rl,EM_max,brite_flag);
else
     % simulate the measurement
    mode_counts = simulateMeasurement(n_pho, p,...
                                             'isPoiss',1,...
                                             'dark_lambda',dark_lambda );
    
    [s_b_trc, s_x_trc, s_y_trc, loglike_trc, count] = EM(mode_counts,num_sources,src_coords,src_brites,prob_fn,X,Y,rl,EM_max,brite_flag,exoplanet_flag);
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
    
    % source colors
    colors = [  [0 0.4470 0.7410]
                [0.8500 0.3250 0.0980]
                [0.9290 0.6940 0.1250]	
                [0.4660 0.6740 0.1880]	
                [0.4940 0.1840 0.5560]
                [0.3010 0.7450 0.9330]	
                [0.6350 0.0780 0.1840]
                ];    
    
    
	% APERTURE
    figs(1) = figure;
    VisualizeAperture(aperture,R_eff);
    
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
            tiledlayout(2,1)
            p_s = prob_fn_measurement(s_x,s_y);
            % make mode counts figure
            nexttile
            stem(1:numel(mode_counts),mode_counts,'filled','k');
            hold on
            for s = 1:num_sources
                mc_s = round(mode_counts .* (p_s(s,:) * s_b(s) ./ sum(s_b.*p_s,1)));
                stem((1:numel(mode_counts)) + .5*(s)/(num_sources), mc_s,'filled','Color',colors(s,:))
                
            end
            hold off
            set(gca,'yscale','log')
            title({'Photon Counting Measurement','Gram-Schmidt Basis',['Total Photons: ',sprintf('%.e',(sum(mode_counts)))]});
            xlabel('mode index')
            ylabel('# photons')
            if exoplanet_flag
                leg = legend(['Total','Star',arrayfun(@(j)['Planet-',num2str(j)],1:num_sources-1,'UniformOutput', 0)]);
            else
                leg = legend(['Total',arrayfun(@(j)['S',num2str(j)],1:num_sources,'UniformOutput', 0)]);    
            end

            title(leg,'Photon Detections by Source')
            n_labels = arrayfun(@num2str,nj,'UniformOutput', 0);
            m_labels = arrayfun(@num2str,mj,'UniformOutput', 0);
            index_labels = strcat(n_labels,repmat({','},[1,num_modes]),m_labels);
            xticks(1:num_modes)
            xticklabels(index_labels)
            
            % make probability distribution figure for each source
            nexttile
            stem(1:num_modes,sum(s_b.*p_s,1),'filled','k')
            hold on
            for s = 1:num_sources
                stem((1:num_modes) + .5*(s)/(num_sources), p_s(s,:), 'filled','Color',colors(s,:))
            end
            hold off
            
            title({'Modal Probability Distributions'});
            xlabel('mode index')
            ylabel('probability')
            if exoplanet_flag
                leg = legend(['Total','Star',arrayfun(@(j)['Planet-',num2str(j)],1:num_sources-1,'UniformOutput', 0)]);
            else
                leg = legend(['Total',arrayfun(@(j)['S',num2str(j)],1:num_sources,'UniformOutput', 0)]);    
            end
            ylim([0,1])
            title(leg,'Mode Detection Probabilities by Source')
            n_labels = arrayfun(@num2str,nj,'UniformOutput', 0);
            m_labels = arrayfun(@num2str,mj,'UniformOutput', 0);
            index_labels = strcat(n_labels,repmat({','},[1,num_modes]),m_labels);
            xticks(1:num_modes)
            xticklabels(index_labels)
            

        case 'Zernike'
            stem(mode_counts);
            set(gca,'yscal','log')
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
    title('Likelihood Convergence','interpreter','latex')
    xlabel('Iteration')
    ylabel('Log Likelihood')
    subplot(2,1,2)
    plot(err_trc./min_sep)
    title('Localization Error Convergence $\epsilon / \Delta_{min}$','interpreter','latex')
    xlabel('Iteration')
    ylabel('Fractional Localization Error')
    

    % ESTIMATE
    figs(5) = figure;
    if exoplanet_flag
        for s = 1:num_sources
            % plot ground truth
            scatter(s_x(s)/rl,s_y(s)/rl,100*log10(s_b(s)./min(s_b) + 1),colors(s,:),'filled'); 
            hold on
        end
        % plot the final constellation estimate
        scatter(s_x_mle(2:end)/rl,s_y_mle(2:end)/rl,100*log10(s_b_mle(2:end)./min(s_b_mle) + 1),'k','filled','square')
        hold off
        names = ['Star',arrayfun(@(j)['Planet-',num2str(j)],1:num_sources-1,'UniformOutput', 0),'Estimate'];
        fig_title = {'Expectation Maximization',basis,...
            ['Brightness Ratio ', sprintf('%0.0e',max(s_b)/min(s_b)),':1'],...
            ['Photons ',sprintf('%0.0e',sum(mode_counts))],...
            ['Minimum Separation $\sigma /',num2str(round(1/min(pdist(scene(:,1:2))))),'$'],...
            ['Localization Error ', sprintf('%0.3g',err/min_sep) ]};
    else
        % plot ground truth
        scatter(s_x/rl,s_y/rl,50*s_b,'filled','black'); 
        hold on;
        % plot the intermediate estimates
        if count > 1
            scatter(s_x_im/rl,s_y_im/rl,50*s_b_im)
        end
        
        % plot the final constellation estimate
        scatter(s_x_mle/rl,s_y_mle/rl,50*s_b_mle,'red','filled','square')
        hold off
        names = cell(1,count+1);
        names(:) = {''}; names(1) = {'Ground Truth'}; names(end) = {'EM Estimate'};
        fig_title = {'Expectation Maximization',basis};
    end
    title(fig_title,'interpreter','latex') 
    xlabel('$x/\sigma$','interpreter','latex');
    ylabel('$y/\sigma$','interpreter','latex');
    xlim([min(X(:))/rl,max(X(:))/rl]); 
    ylim([min(Y(:))/rl,max(Y(:))/rl]); 
    xticks(linspace(min(X(:))/rl,max(X(:))/rl,5));
    yticks(linspace(min(Y(:))/rl,max(Y(:))/rl,5));
    legend(names)
    axis square
    
end

end