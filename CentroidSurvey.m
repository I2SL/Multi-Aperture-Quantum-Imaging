function CentroidSurvey(array_id,num_workers)

    addpath('utils/')

    % get the job array id
    if ischar(array_id)
        array_id = str2double(array_id);
    end
    
    % set the random number generator seed
    rng(array_id)
    
    % make the DS structure
    DS = DSformat();
    
    % make the save directory
    mkdir(DS.save_dir)
    
    % get configuration indices
    [a,n,m,p,b] = ind2sub(DS.cfg_size,array_id);
    cfg_id = {a,n,m,p,b};
    
    % parfor configuration variables
    trials = DS.trials;
    EM_iters_max = DS.EM_iters_max;         % number of EM iterations to use per initialization
    EM_cycles = DS.EM_cycles;               % number of EM initializations to run per measurement
    aperture = DS.apertures{a};
    num_src = DS.num_src(n);
    min_sep_frac = DS.min_sep_frac(m);
    num_pho = DS.num_pho(p);
    basis = DS.basis{b};
    max_order = DS.max_order;
   
    % multi-aperture parameters
    ap_num = size(aperture,1);
    aper_coords = aperture(:,1:2);
    aper_rads = aperture(:,3); 

    % get the effective aperture diameter
    if ap_num>1
        B = squareform(pdist(aper_coords));             % sub-aperture pairwise centroid distances
        D = aper_rads + aper_rads';                     % sub-aperture pairwise radius sums
        assert(all(triu(B,1) >= triu(D,1),'all'));      % check if any apertures overlap

        % set the effective aperture diameter to the minimum enclosing circle diameter
        cm_coords = aper_coords - mean(aper_coords,1);                          % get centered coordinates (with respect to centroid of centroids -- still need to prove with this works)
        tangent_pts = cm_coords + cm_coords.*aper_rads./vecnorm(cm_coords,2,2); % candidate tangent points where the enclosing circle might touch
        [~,~,R_eff] = SmallestEnclosingCircle(tangent_pts(:,1)',tangent_pts(:,2)'); % effective aperture radius
        D_eff = 2*R_eff;
    else
        R_eff = aper_rads(1);
        D_eff = 2*R_eff;                         % set the effective aperture diameter to that of the input aperture
    end


    % The rayleigh length of the multi-aperture system is defined to be the
    % same rayleigh length as a single circular hard aperture with diameter
    % D_eff, where D_eff is the diameter of the minimum enclosing circle 
    % that contains all of the sub-apertures.
    rl = 2*pi * 1.2197/D_eff; % rayleigh length in units of [rads/length] 
    
    % image plane discretization
    [X,Y] = meshgrid(rl * linspace(-.5,.5,DS.pos_samp));       % image plane coordinates [rad/u]

    % make sure min source separation is greater than the resolution of the image plane
    dx = X(1,2) - X(1,1);
    if min_sep_frac < dx/rl
        warning('Image plane discretization is coarser than minimum source separation')
    end
    
    
    % setup basis
    switch basis
        case 'Gram-Schmidt'
            
            % number of modes
            num_modes = (max_order+1)^2;

            % aperture plane discretization
            [Kx,Ky,d2k] = ApertureKxKy(aper_coords,aper_rads,DS.mom_samp);     % Kx [u], Ky [u], d2k [u^2]
            
            % assign A_tot to be the area of tbe discretized aperture
            A_tot = numel(Kx)*d2k;            

            % Create Gram-Schmidt basis
            GS_basis_mom = genGramSchmidtBasis_mom(DS.max_order,Kx,Ky,d2k);              % basis functions in momentum space
            GS_basis_pos = Basis_GramSchmidt_mom([X(:),Y(:)],Kx,Ky,d2k,GS_basis_mom);    % basis functions in position space

            % probability function
            GS_basis_pos = reshape(GS_basis_pos,[size(X),size(GS_basis_pos,2)]);                    % 2D basis matrix stack
            prob_fn = @(xq,yq) ModalProb_GramSchmidt_pos([xq,yq],X,Y,GS_basis_pos,A_tot);


        case 'Zernike'
            
            % number of modes
            num_modes = ap_num * (max_order+1)*(max_order+2)/2;

            % indices
            [nj,mj,vj] = Indices_MixedAperture(max_order,ap_num);

            % Create Mixed-Aperture Zernike Basis
            U = dftmtx(ap_num)/sqrt(ap_num);   % a unitary matrix for mixing aperture-local modes

            % probability function handle
            prob_fn = @(xq,yq) ModalProb_MixedAperture([xq,yq],nj,mj,vj,U,aperture);

        case 'Direct-Detection'
            
            % number of modes
            num_modes = numel(X); 

            % probability function
            prob_fn = @(xq,yq) ModalProb_DirectImaging([xq,yq],X,Y,aperture);

    end
    
    
    disp(['-------Configuration: ' num2str(array_id),'/',num2str(prod(DS.cfg_size)),'--------'])    
    
    % run parameter scans using matlab's Parallel Computing Toolbox
    parpool(num_workers)
    
    %for t=1:DS.trials
    parfor t=1:DS.trials
        
        % configure centroid-aligned scene
        align_centroid = 1;
        s_b = ones(num_src, 1) / num_src;    % relative source brightnesses
        src_coords_frac = genMinDistConstellation(s_b,min_sep_frac,align_centroid); % source coordinates in fractional rayleigh units [rl]
        aligned_scene = [src_coords_frac, s_b];
               
        % model the centroid misalignment as a gaussian random variable with
        % truncated support. PDF is 0 beyond the rayleigh FOV.
        align_sigma_frac = 0.1; % standard deviation of misalignment [rl] 
        centroid_shift = mvnrnd([0,0],align_sigma_frac*eye(2));
        while norm(centroid_shift) > 0.5
            centroid_shift = mvnrnd([0,0],align_sigma_frac*eye(2));
        end
        
        % the unaligned scene
        unaligned_scene = [src_coords_frac + centroid_shift, s_b];
                
        % source coordinates of unaligned scene in position-space units
        unaligned_src_coords = rl*unaligned_scene(:,1:2);            % source coordinates [rad/length]
              
        % sample the number of collected photons
        N = poissrnd(num_pho);                  % total number of collected photons
        N1 = round(N/2);                        % photons allocated for centroid pre-estimate (Grace et. al. 2020 limit for deep sub-rayleigh)
        
        % direct-imaging probability over unaligned source distribution
        lam = D_eff/min(aper_rads);     % only consider PSF up to the first zero of the largest bessel function (the rayleigh limit of the smallest aperture)
        [X_DD,Y_DD] = meshgrid(lam*rl*linspace(-.5,.5,round(lam*DS.pos_samp))); % make the domain over which we perform direct imaging
        p_DD = sum(s_b .* ModalProb_DirectImaging(unaligned_src_coords,X_DD,Y_DD,aperture),1);
        
        % check if the mode and mean of the PSF are the same. Otherwise the
        % mean photon location is not the MLE
        [~,mo_id] = max(p_DD);
        mo_psf = [X_DD(mo_id),Y_DD(mo_id)];
        mu_psf = p_DD/sum(p_DD)*[X_DD(:),Y_DD(:)];
        tolerance = 0.01;
        if norm(mo_psf - mu_psf) > tolerance
            warning('Centroid estimator may not be the MLE. PSF mode and mean are unequal')
        end
        
        % collect N1 Direct Imaging photons from the unaligned scene to
        % estimate the centroid
        [pho_xy_id,~] = simulateMeasurement(N1,p_DD,0);
        x_DD = X_DD(pho_xy_id);
        y_DD = Y_DD(pho_xy_id);
        
        %pho_xy_id = datasample(1:numel(p_DD),N1,'Weights',p_DD);  % get the cell indices in which each photon was detected
        centroid_est = mean([x_DD',y_DD']/rl,1);      % estimate the centroid [in rayleigh units]
        
        % correct the unaligned scene by recentering to the centroid estimate
        corrected_scene = [(centroid_shift-centroid_est) + src_coords_frac,s_b];
        
        % if any of the sources in the corrected scene fall outside the
        % FOV, then
        % outside of the rayleigh FOV 
        if any(vecnorm(corrected_scene(:,1:2),2,2) > 0.5)
            u = (centroid_shift-centroid_est) / norm(centroid_shift-centroid_est);
            delta = (.5-max(vecnorm(src_coords_frac,2,2))) * (1-1e-4)*u ;
            corrected_scene  = [delta+src_coords_frac,s_b];
        end
        
        assert( all(vecnorm(corrected_scene(:,1:2),2,2) < 0.5))
        
        % group the misaligned, corrected, and aligned scenes
        scene_group = zeros(num_src,3,1,3);
        scene_group(:,:,1,1) = unaligned_scene;
        scene_group(:,:,1,2) = corrected_scene;
        scene_group(:,:,1,3) = aligned_scene;
        
        %{
        % visualize
        figure
        hold on
        colormap gray
        d2x_DD = (X_DD(1,2)-X_DD(1,1)).^2;
        imagesc(lam*[-.5,.5],lam*[.5,-.5],flipud(reshape(p_DD/d2x_DD,size(X_DD))))
        cbar = colorbar;
        scatter(X_DD(pho_xy_id)'/rl,Y_DD(pho_xy_id)'/rl,10,'filled','yellow','MarkerEdgeColor','black')
        scatter(centroid_shift(1),centroid_shift(2),36,'+blue');
        scatter(centroid_est(1),centroid_est(2),36,'+red');
        hold off
        axis equal
        xlabel('x [rl]')
        ylabel('y [rl]')
        xlim(lam*[-.5,.5])
        ylim(lam*[-.5,.5])
        legend({'Photon Arrivals','Centroid','Centroid Estimate'})
        ylabel(cbar,'Photon Incidence Probability Density') 
        title({'Centroid Estimation with Direct Detection',['Direct Detection Photons: ',num2str(N1)]})
         
        % show the different scenes
        figure
        hold on
        scatter(unaligned_scene(:,1),unaligned_scene(:,2),'filled','blue')
        scatter(corrected_scene(:,1),corrected_scene(:,2),'filled','red')
        scatter(aligned_scene(:,1),aligned_scene(:,2),'filled','black')
        hold off
        xlim([-.5,.5])
        ylim([-.5,.5])
        xlabel('x [rl]')
        ylabel('y [rl]')
        title('Target Scenes')
        legend({'Misaligned Target','Centroid Pre-Estimate Correction','Perfectly Aligned'})
        %}
        
        % a container to hold measurements for each scene's centroid config
        measurement_group = zeros(1,num_modes,1,3);
        
        % perform modal imaging on all scenes in the scene group
        for s = 1:3
            
            % get the number of photons allocated for modal imaging
            if s==2
                N2 = N-N1;
            elseif s==1 || s==3
                N2 = N;
            end
            
            % get the scene
            scene = scene_group(:,:,s);
            s_x = rl* scene(:,1);
            s_y = rl* scene(:,2);
            src_coords = [s_x,s_y];
            
            % get modal probabilities for the given source distribution
            p_scene = sum(s_b .* prob_fn(s_x,s_y),1);
            
             % simulate the modal measurement
            [~, mode_counts] = simulateMeasurement(N2, p_scene, 0);
            measurement_group(1,:,1,s) = mode_counts;
            
            % run expectation maximization on the measurement
            est_scene = zeros(num_src,3,EM_cycles,3);
            err = zeros(1,1,EM_cycles,3);
            loglike = zeros(1,1,EM_cycles,3);
            EM_iters = zeros(1,1,EM_cycles,3);
            
            % different EM initializations
            tic
            for k = 1:EM_cycles 
                
                % fine MLE of scene parameters given the measurement
                [s_b_trc, s_x_trc, s_y_trc, loglike_trc, iters] = EM(mode_counts,num_src,prob_fn,X,Y,rl,EM_iters_max,0);
                
                % final scene parameter estimates
                s_b_mle = s_b_trc(:,end); s_x_mle = s_x_trc(:,end); s_y_mle = s_y_trc(:,end);
                
                % mle scene estimate
                est_coords = [s_x_mle,s_y_mle];
                est_coords_frac = est_coords/rl;
                est_brites = s_b_mle;
                
                % collect results of EM cycle
                est_scene(:,:,k,s) = [est_coords_frac, est_brites];   % scene estimate
                loglike(1,1,k,s) = loglike_trc(1,end);                    % log likelihood
                err(1,1,k,s) = LocalizationError(src_coords, est_coords); % localization error
                EM_iters(1,1,k,s) = iters;                                % EM iterations
            end 
            toc
            disp([num2str(DS.EM_cycles),' EM cycles completed.'])

        end
        
        % data stucture for trial
        cfg_data(t).rl = rl;
        cfg_data(t).N1 = N1;
        cfg_data(t).N2 = N2;
        cfg_data(t).DD_centroid_pho_xy = [x_DD',y_DD'];
        cfg_data(t).centroid = centroid_shift;
        cfg_data(t).centroid_est = centroid_est;
        cfg_data(t).measurement = measurement_group;
        cfg_data(t).scene = scene_group;
        cfg_data(t).scene_est = est_scene;
        cfg_data(t).loglike = loglike;
        cfg_data(t).err = err;
        cfg_data(t).EM_iters = EM_iters;  

        % display trial completion
        disp(['Trials Completed: ',num2str(t),'/',num2str(trials)])
    end
    
    DS.data(cfg_id{:}) = {cfg_data};

    % save current data structure
    fname = [num2str(array_id),'_cfg','.mat'];
    save(fullfile(DS.save_dir,fname),'cfg_id','DS')
    
end
    