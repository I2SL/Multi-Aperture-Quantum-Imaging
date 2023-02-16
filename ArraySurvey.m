function ArraySurvey(array_id)

    addpath('utils/')

    % get the job array id
    if ischar(array_id)
        array_id = str2double(array_id);
    end
    
    % make the DS structure
    DS = DSformat();
    
    % get configuration indices
    [b,p,a,n,m] = ind2sub(DS.cfg_size,array_id);
    cfg_id = {b,p,a,n,m};
    
    % parfor configuration variables
    EM_iters = DS.EM_iters;             % number of EM iterations to use per initialization
    EM_cycles = DS.EM_cycles;           % number of EM initializations to run per measurement
    basis = DS.basis{b};
    num_pho = DS.num_pho(p);
    aper_coords = DS.apertures{a};
    num_src = DS.num_src(n);
    min_sep_frac = DS.min_sep_frac(m);
                       
    % rescale aperture coordinates to be in the reference unit
    aper_coords = aper_coords/DS.ref_unit;     % sub-aperture coordinates [u]
        
    % multi-aperture parameters
    ap_num = size(aper_coords,1);           % number of sub-apertures
    A_sub = pi*DS.subap_radius^2;              % subaperture collection area [u^2]
    A_tot = ap_num * A_sub;                 % total collection area of the multi-aperture system [u^2]
    if ap_num > 1
        B = pdist(aper_coords);                 % baseline lengths [u]
        max_B = max([1,B]);                     % max baseline [u]
    else
        max_B = 1;
    end

    % rayleigh lengths
    rl_sub = 2*pi*1.22;                     % sub-aperture rayleigh length [rad/u]
    rl = rl_sub/max_B;                      % effective aperture rayleigh length [rad/u]

    % image plane discretization
    [X,Y] = meshgrid(rl * linspace(-.5,.5,DS.img_samp));       % image plane coordinates [rad/u]
    dx = X(1,2) - X(1,1);
    
    % setup probability functions for the basis
    switch basis
        case 'Gram-Schmidt'

            % indices
            [nj,~] = Indices_GramSchmidt(DS.max_order);

            % number of modes
            num_modes = numel(nj);

            % aperture plane discretization
            [Kx,Ky,d2k] = ApertureKxKy(aper_coords,DS.subap_samp);     % Kx [u], Ky [u], d2k [u^2]

            % Create Gram-Schmidt basis
            GS_basis_mom = genGramSchmidtBasis_mom(DS.max_order,Kx,Ky,d2k);                 % basis functions in momentum space
            GS_basis_pos = Basis_GramSchmidt_mom([X(:),Y(:)],Kx,Ky,d2k,GS_basis_mom);    % basis functions in position space
            GS_basis_pos = reshape(GS_basis_pos,[size(X),num_modes]);                    % 2D basis matrix stack

            % probability function
            prob_fn = @(xq,yq) ModalProb_GramSchmidt_pos([xq,yq],X,Y,GS_basis_pos,A_tot);


        case 'Zernike'

            % indices
            [nj,mj,vj] = Indices_MixedAperture(DS.max_order,ap_num);

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
    
    
    %%
    
    % run parameter scans using matlab's Parallel Computing Toolbox
    parpool(94)
    
    disp(['-------Configuration: ' num2str(array_id),'/',num2str(prod(DS.cfg_size)),'--------'])    
    
    %for t=1:DS.trials
    parfor t=1:DS.trials
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % configure scene
        src_coords_frac = genMinDistConstellation(num_src,min_sep_frac); % source coordinates in fractional rayleigh units [x/rl]
        src_brites = ones(num_src,1) / num_src;
        scene = [src_coords_frac, src_brites];

        % source distribution
        src_coords = rl*scene(:,1:2);            % source coordinates [rad/u]
        s_x = src_coords(:,1); s_y = src_coords(:,2); 
        s_b = scene(:,3);                        % relative source brightnesses

        % make sure min source separation is greater than the resolution of the image plane
        if rl*min_sep_frac < dx
            warning('Image plane discretization is coarser than minimum source separation')
        end

        % get modal probabilities for the given source distribution
        p_scene = sum(s_b .* prob_fn(s_x,s_y),1);
        
        % simulate the measurement
        [~, mode_counts] = simulateMeasurement(num_pho, p_scene);
        
        
        % different EM initializations
        est_scene = zeros(num_src,3,EM_cycles);
        err = zeros(1,1,EM_cycles);
        loglike = zeros(1,1,EM_cycles);
        for k = 1:EM_cycles   
            tic         
            % find MLE of scene parameters given the measurement
            [s_b_trc, s_x_trc, s_y_trc, loglike_trc,~] = EM(mode_counts,num_src,prob_fn,X,Y,rl,EM_iters);

            % final scene parameter estimates
            s_b_mle = s_b_trc(:,end); s_x_mle = s_x_trc(:,end); s_y_mle = s_y_trc(:,end);
            
            % mle scene estimate
            est_coords = [s_x_mle,s_y_mle];
            est_coords_frac = est_coords/rl;
            est_brites = s_b_mle;

            % collect results of EM cycle
            est_scene(:,:,k) = [est_coords_frac, est_brites];   % scene estimate
            loglike(k) = loglike_trc(1,end);                    % log likelihood
            err(k) = LocalizationError(src_coords, est_coords); % localization error
            toc
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % add to data element
        cfg_data(t,:) = {rl,scene,mode_counts,est_scene,loglike,err};    

        % display trial completion
        disp(['Trials Completed: ',num2str(t),'/',num2str(DS.trials)])
    end
    
    DS.data(cfg_id{:}) = {cfg_data};

    % save current data structure
    fname = [num2str(array_id),'cfg_Survey_2-7-23','.mat'];
    save(fullfile('data_out',fname),'cfg_id','DS')    
end