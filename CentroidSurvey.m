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
    [b,p,a,n,m] = ind2sub(DS.cfg_size,array_id);
    cfg_id = {b,p,a,n,m};
    
    % parfor configuration variables
    trials = DS.trials;
    EM_iters_max = DS.EM_iters_max;             % number of EM iterations to use per initialization
    EM_cycles = DS.EM_cycles;               % number of EM initializations to run per measurement
    basis = DS.basis{b};
    num_pho = DS.num_pho(p);
    aperture = DS.apertures{a};
    num_src = DS.num_src(n);
    min_sep_frac = DS.min_sep_frac(m);
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

            % assign A_tot to be the true aperture area
            A_tot = sum(pi*aper_rads.^2);

            % indices
            [nj,mj,vj] = Indices_MixedAperture(max_order,ap_num);

            % Create Mixed-Aperture Zernike Basis
            U = dftmtx(ap_num)/sqrt(ap_num);   % a unitary matrix for mixing aperture-local modes

            % probability function handle
            prob_fn = @(xq,yq) ModalProb_MixedAperture([xq,yq],nj,mj,vj,U,aperture);

        case 'Direct-Detection'

            % probability function
            prob_fn = @(xq,yq) ModalProb_DirectImaging([xq,yq],X,Y,aperture);

    end
    
    
    
    % run parameter scans using matlab's Parallel Computing Toolbox
    %parpool(num_workers)
    
    disp(['-------Configuration: ' num2str(array_id),'/',num2str(prod(DS.cfg_size)),'--------'])    
    
    for t=1:DS.trials
    %parfor t=1:DS.trials
        
        % configure centroid-aligned scene
        align_centroid = 1;
        s_b = ones(num_src, 1) / num_src;    % relative source brightnesses
        src_coords_frac = genMinDistConstellation(s_b,min_sep_frac,align_centroid); % source coordinates in fractional rayleigh units [rl]
        aligned_scene = [src_coords_frac, s_b];
        
        % configure centroid-unaligned scene with a random displacement
        % that stays within the FOV
        delta_r_max = max(0.5-vecnorm(src_coords_frac,2,2));
        [c_x,c_y] = pol2cart(2*pi*rand(1),delta_r_max*rand(1));
        centroid_shift = [c_x,c_y]; 
        unaligned_scene = [src_coords_frac + centroid_shift, s_b];
                
        % source positions and brightneses
        unaligned_src_coords = rl*unaligned_scene(:,1:2);            % source coordinates [rad/length]
              
       % sample the number of collected photons
        N = poissrnd(num_pho);                  % total number of collected photons
        N1 = round(sqrt(N));                    % photons allocated for centroid pre-estimate
        
        % direct-imaging probability over unaligned source distribution
        %[X_DD,Y_DD] = meshgrid(D_eff*rl*linspace(-.5,.5,D_eff*DS.pos_samp)); % make the domain over which we perform direct imaging
        X_DD = X; Y_DD = Y;
        p_DD = sum(s_b .* ModalProb_DirectImaging(unaligned_src_coords,X_DD,Y_DD,aperture),1);
        
        % collect N1 Direct Imaging photons from the unaligned scene to
        % estimate the centroid
        [pho_xy_id,~] = simulateMeasurement(N1,p_DD,0);
        
        %pho_xy_id = datasample(1:numel(p_DD),N1,'Weights',p_DD);  % get the cell indices in which each photon was detected
        centroid_est = mean([X_DD(pho_xy_id)',Y_DD(pho_xy_id)']/rl,1);      % estimate the centroid [in rayleigh units]
        
        % if centroid estimate pushes sources outside of rayleigh FOV then
        % rescale it back into the FOV
        if norm(centroid_est)>delta_r_max
            centroid_est = centroid_est / norm(centroid_est) * delta_r_max;
        end
        
        %{
        % visualize
        figure
        hold on
        colormap gray
        d2x_DD = (X_DD(1,2)-X_DD(1,1)).^2;
        imagesc(D_eff*[-.5,.5],D_eff*[.5,-.5],flipud(reshape(p_DD/d2x_DD,size(X_DD))))
        cbar = colorbar;
        scatter(X_DD(pho_xy_id)'/rl,Y_DD(pho_xy_id)'/rl,10,'filled','yellow')
        scatter(centroid_shift(1),centroid_shift(2),36,'+blue');
        scatter(centroid_est(1),centroid_est(2),36,'+red');
        hold off
        axis equal
        xlabel('x [rl]')
        ylabel('y [rl]')
        xlim(D_eff*[-.5,.5])
        ylim(D_eff*[-.5,.5])
        legend({'Photon Arrivals','Centroid','Centroid Estimate'})
        ylabel(cbar,'Photon Incidence Probability Density') 
        title({'Centroid Estimation with Direct Detection',['Direct Detection Photons: ',num2str(N1)]})
        %}
        
        % correct the unaligned scene by recentering to the centroid estimate
        corrected_scene = [src_coords_frac+centroid_shift-centroid_est,s_b];
        
        % group the misaligned, corrected, and aligned scenes
        scene_group = zeros(num_src,3,3);
        scene_group(:,:,1) = unaligned_scene;
        scene_group(:,:,2) = corrected_scene;
        scene_group(:,:,3) = aligned_scene;
        
        
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
            
            % run expectation maximization on the measurement
            est_scene = zeros(num_src,3,EM_cycles);
            err = zeros(1,1,EM_cycles);
            loglike = zeros(1,1,EM_cycles);
            EM_iters = zeros(1,1,EM_cycles);
            
            % different EM initializations
            for k = 1:EM_cycles 
                
                % fine MLE of scene parameters given the measurement
                [s_b_trc, s_x_trc, s_y_trc, loglike_trc, iters] = EM(mode_counts,num_src,prob_fn,X,Y,rl,EM_iters_max,0);
                
                % final scene parameter estimates
                s_b_mle = s_b_trc(:,end); s_x_mle = s_x_trc(:,end); s_y_mle = s_y_trc(:,end);
                
                % mle scene estimate
                est_coords = [s_x_mle,s_y_mle];
                est_coords_frac = est_coords/rl;
                est_brites = s_b_mle;
                
                % FIGURE OUT HOW WE WILL BE STORING THE RESULTS
                % collect results of EM cycle
                est_scene(:,:,k) = [est_coords_frac, est_brites];   % scene estimate
                loglike(k) = loglike_trc(1,end);                    % log likelihood
                err(k) = LocalizationError(src_coords, est_coords); % localization error
                EM_iters(k) = iters;                                % EM iterations
            end           
        end
        disp([num2str(DS.EM_cycles),' EM cycles completed.'])
        toc
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % add to data element
        cfg_data(t,:) = {rl,scene,mode_counts,est_scene,loglike,err,EM_iters};    

        % display trial completion
        disp(['Trials Completed: ',num2str(t),'/',num2str(trials)])
    end
    
    DS.data(cfg_id{:}) = {cfg_data};

    % save current data structure
    fname = [num2str(array_id),'cfg_2-23-23','.mat'];
    save(fullfile(DS.save_dir,fname),'cfg_id','DS')
    
end
    