function DS = DSformat()
    % DS.cfg_data = {rl, scene, measurement, mle_scene, likelihood, err, EM_iterations};
    % cfg = [b,p,a,n,m,t]   --> configuration index vector [b = basis index, p = num_pho index, a = aperture index, n = num_src index, m = min_sep index, t = cfg_trials index]
    % rl                    --> the rayleigh length of the configuration
    % scene                 --> [n_src X 3 ] array contianing the scene parameters : x = scene(:,1), y = scene(:,2), b = scene(:,3) 
    % measurement           --> [1 X n_modes ] array of photon counts in each mode
    % mle_scene             --> [n_src x 3 x EM_cycles] array maximum likelihood estimates arrived at for the EM algorithm 
    % likelihood            --> [1 x 1 x EM_cycles] array with the likelihoods of the estimate produced at each EM cycle
    % err                   --> [1 x 1 x EM_cycles] array with the fractional localization error of the estimate produced at each EM cycle

    % save directory
    save_dir = fullfile('Survey_2-27-23_Centroid_Dependence','data_out');
    
    
    % constants
    trials = 50;        % trials per configuration
    mom_samp = 67;      % sample density of aperture plane [Gram-Schmidt] [samples/length]
    pos_samp = 129;     % image plane samples (must be odd!)
    EM_iters_max = 100; % max EM iterations
    EM_cycles = 15;     % number of times to run the EM algorithm (with different initializations) on the same measurement
    max_order = 5;      % max basis order for GS and Zernike    
    
    % setup apertures
    D_eff = 30;         % multi-aperture effective diameter [length]
    R_eff = D_eff/2;    % multi-aperture effective radius   [length]
    d = 3;              % sub-aperture diameter             [length]
    r = d/2;            % sub-apeture radius                [length]
    ap3 = [Polygon(3,0,'radius',R_eff-r),r*ones(3,1)];
    ap9 = [Polygon(9,0,'radius',R_eff-r),r*ones(9,1)];
    golay9 = [Golay9(R_eff-r),r*ones(9,1)];    
    apertures = {ap3,ap9,golay9};
    aperture_names = {'3 Aperture','9 Aperture','Golay-9'};    

    % data structure with some limited functionality
    DS = struct();
    
    % const properties
    DS.timestamp = datetime;
    DS.save_dir = save_dir;                 % save directory for the data structure
    DS.trials = trials;                     % how many unique scenes/measurements we generate per configuration
    DS.EM_iters_max = EM_iters_max;         % max number of iterations EM is allowed to run for
    DS.EM_cycles = EM_cycles;               % how many times we instantiate EM on the same measurement
    DS.pos_samp = pos_samp;                 % image plane samples (should be odd)
    DS.mom_samp = mom_samp;                 % sub aperture samples (should be odd) [for gram-schmidt basis mainly]
    DS.R_eff = R_eff;                       % effective aperture radius
    DS.max_order = max_order;               % max modal order 
    DS.cfg_idx_names = {'Basis Index (b)','Photon Number Index (p)','Aperture Index (a)','Source Number Index (n)','Min Separation Index (m)'};
    DS.cfg_data_names = {'Rayleigh Length','Scene','Measurement Mode Counts','Estimated Scene','Log Likelihood','Error'};
    
    % array properties (parameter scans)
    DS.num_pho = [2e3,2e4,2e5];
    DS.basis = {'Gram-Schmidt','Direct-Detection'};
    DS.min_sep_frac = 2.^(linspace(-6,-3,7)); % fractional rayleigh units
    DS.apertures = apertures;                 
    DS.aperture_names = aperture_names;
    DS.num_src = 3:5;
    DS.cfg_size = [numel(DS.basis),numel(DS.num_pho),numel(DS.apertures),numel(DS.num_src),numel(DS.min_sep_frac)];  % the dimensionality of the parameter space range
    DS.data = cell(DS.cfg_size);
    
end


