function DarkCurrentSurvey(array_id,num_workers)

    addpath('utils/')

    % get the job array id
    if ischar(array_id)
        array_id = str2double(array_id);
    end
    
    % make the DS structure
    DS = DSformat_dark_current();
    
    % make the save directory
    mkdir(DS.save_dir)
    
    % get configuration indices
    [a,n,m,p,b,e1,e2] = ind2sub(DS.cfg_size,array_id);
    cfg_id = {a,n,m,p,b,e1,e2};
    
    % parfor configuration variables
    trials = DS.trials;
    EM_iters_max = DS.EM_iters_max;         % number of EM iterations to use per initialization
    EM_cycles = DS.EM_cycles;               % number of EM initializations to run per measurement
    aperture = DS.apertures{a};
    rl = DS.rl(a);                          % rayleigh length in units of [rads/length] 
    num_src = DS.num_src(n);
    min_sep_frac = DS.min_sep_frac(m);
    num_pho = DS.num_pho(p);
    basis = DS.basis{b};
    max_order = DS.max_order;
    dark_lambda = DS.dark_lambda(e1);
    phase_sigma = DS.phase_sigma(e2);
    
    % change the number of EM cycles to 1 for DD. All initializations
    % converge to the same result so no need to do multiple trials on the
    % same measurement.
    if strcmp(basis,'Direct-Detection')
        EM_cycles = 1;
        if DS.dark_lambda(e1) ~= 0
            return
        end
        
        if DS.phase_sigma(e2) ~= 0
            return
        end
    end
    
    % monolith cannot have relative phasing error
    if strcmp(DS.aperture_names{a},'Monolith')
        if DS.phase_sigma(e2) ~= 0
            return
        end
    end
   
   
    % multi-aperture parameters
    ap_num = size(aperture,1);
    aper_coords = aperture(:,1:2);
    aper_rads = aperture(:,3); 
    D_eff = 2*pi * 1.2197 / rl; % Effective aperture diameter [length]
    
    % Gram-Schmidt basis image pland discretization (only has to be defined
    % over sub-rayleigh regime if the scene is sub-rayleigh)
    [X_GS,Y_GS] = meshgrid(rl * linspace(-.5,.5,DS.pos_samp));       % image plane coordinates [rad/u]
    
    % direct detection image plane discretization
    samp_scale = D_eff/(2*min(aper_rads));     % only consider PSF up to the first zero of the largest bessel function (the rayleigh limit of the smallest aperture) which covers most of the probability in the PSF
    [X_DD,Y_DD] = meshgrid(samp_scale*rl*linspace(-.5,.5,round(samp_scale*DS.pos_samp))); % make the domain over which we perform direct imaging

    % make sure min source separation is greater than the resolution of the image plane
    dx_GS = X_GS(1,2) - X_GS(1,1);
    dx_DD = X_DD(1,2) - X_DD(1,1);
    if min_sep_frac < dx_GS/rl || min_sep_frac < dx_DD/rl
        warning('Image plane discretization is coarser than minimum source separion')
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
            GS_basis_pos = Basis_GramSchmidt_mom([X_GS(:),Y_GS(:)],Kx,Ky,d2k,GS_basis_mom);    % basis functions in position space

            % probability function
            GS_basis_pos = reshape(GS_basis_pos,[size(X_GS),size(GS_basis_pos,2)]);                    % 2D basis matrix stack
            prob_fn = @(xq,yq) ModalProb_GramSchmidt_pos([xq,yq],X_GS,Y_GS,GS_basis_pos,A_tot);

            prob_fn_measurement = prob_fn;

        case 'Zernike'
            
            % number of modes
            num_modes = ap_num * (max_order+1)*(max_order+2)/2;

            % indices
            [nj,mj,vj] = Indices_MixedAperture(max_order,ap_num);

            % Create Mixed-Aperture Zernike Basis
            U = dftmtx(ap_num)/sqrt(ap_num);   % a unitary matrix for mixing aperture-local modes

            % probability function handle
            prob_fn = @(xq,yq) ModalProb_MixedAperture([xq,yq],nj,mj,vj,U,aperture);
            
            U_noisy = U+exp(1i*2*pi*(normrnd(0,phase_sigma,size(U))));
            prob_fn_measurement = @(xq,yq) ModalProb_MixedAperture([xq,yq],nj,mj,vj,U_noisy,aperture);
            
        case 'Direct-Detection'
            
            % number of modes
            num_modes = numel(X_DD); 

            % probability function
            prob_fn = @(xq,yq) ModalProb_DirectImaging([xq,yq],X_DD,Y_DD,aperture);
            prob_fn_measurement = prob_fn;
    end
    
    
    disp(['-------Configuration: ' num2str(array_id),'/',num2str(prod(DS.cfg_size)),'--------'])    
    
    
    %for t=1:DS.trials
    
    % run parameter scans using matlab's Parallel Computing Toolbox
    parpool(num_workers)
    parfor t=1:trials
        
        % set the random number generator seed (must include the parallel
        % worker in the seed for different instances to develop different scenes)
        rng(array_id*(trials-1) + t)
    
        % configure centroid-aligned scene
        align_centroid = 1;
        s_b = ones(num_src, 1) / num_src;    % relative source brightnesses
        src_coords_frac = genMinDistConstellation(s_b,min_sep_frac,align_centroid); % source coordinates in fractional rayleigh units [rl]
        scene = [src_coords_frac, s_b];     % scene object      
        
        % run expectation maximization on the measurement
        est_scene = zeros(num_src,3,EM_cycles);
        err = zeros(1,1,EM_cycles);
        loglike = zeros(1,1,EM_cycles);
        EM_iters = zeros(1,1,EM_cycles);        
                   
        % get the scene
        s_x = rl * scene(:,1);
        s_y = rl * scene(:,2);
        src_coords = [s_x,s_y];
        
        % mean photon number
        N = num_pho;
        
        % simulate the measurement
        switch basis
            case 'Direct-Detection'
                % sample photon arrivals from the probability
                % distribution given by the ground truth sources and
                % the aperture configuration
                p_DD = sum(s_b .*prob_fn(s_x,s_y),1);
                [pho_xy_id, mode_counts] = simulateMeasurement(N, p_DD,'isPoiss',1);                 
                x_DD = X_DD(pho_xy_id);
                y_DD = Y_DD(pho_xy_id);
            otherwise
                % get modal probabilities for the given source distribution
                p_modes = sum(s_b .* prob_fn(s_x,s_y),1);
                [~, mode_counts] = simulateMeasurement(N, p_modes, 'isPoiss',0,'dark_lambda',dark_lambda);   
        end

        % different EM initializations
        tic
        for k = 1:EM_cycles 

            % find MLE of scene parameters given the measurement
            switch basis 
                case 'Direct-Detection'
                    [s_b_trc, s_x_trc, s_y_trc, loglike_trc, iters] = EM_DD([x_DD',y_DD'],num_src,aperture,rl,EM_iters_max,0);                        
                otherwise
                    [s_b_trc, s_x_trc, s_y_trc, loglike_trc, iters] = EM(mode_counts,num_src,src_coords,prob_fn,X_GS,Y_GS,rl,EM_iters_max,0);                        
            end

            % final scene parameter estimates
            s_b_mle = s_b_trc(:,end); 
            s_x_mle = s_x_trc(:,end); 
            s_y_mle = s_y_trc(:,end);

            % mle scene estimate
            est_coords = [s_x_mle,s_y_mle];
            est_coords_frac = est_coords/rl;
           
            
            % collect results of EM cycle
            est_scene(:,:,k) = [est_coords_frac, s_b_mle];                              % scene estimate
            loglike(1,1,k) = loglike_trc(1,end);                                        % log likelihood
            err(1,1,k) = LocalizationError(src_coords, est_coords)/(rl*min_sep_frac);   % average fractional localization error
            EM_iters(1,1,k) = iters;                                                    % EM iterations
        end 
        toc
        disp([num2str(EM_cycles),' EM cycles completed.'])

        
        % data stucture for trial
        cfg_data(t).N = N;
        cfg_data(t).scene = scene;
        cfg_data(t).est_scene = est_scene;
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
    