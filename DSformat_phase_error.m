function DS = DSformat_phase_error()
    % DS.cfg_data = {rl, scene, measurement, mle_scene, likelihood, err, EM_iterations};
    % cfg = [a,n,m,p,b] --> configuration index vector [a = aperture index, n = num_src index, m = min_sep index, p = num_pho index, b = basis index]
    %
    % rl                    --> the rayleigh length of the configuration
    % scene                 --> [n_src X 3 ] array contianing the scene parameters : x = scene(:,1), y = scene(:,2), b = scene(:,3) 
    % measurement           --> [1 X n_modes ] array of photon counts in each mode
    % mle_scene             --> [n_src x 3 x EM_cycles] array maximum likelihood estimates arrived at for the EM algorithm 
    % likelihood            --> [1 x 1 x EM_cycles] array with the likelihoods of the estimate produced at each EM cycle
    % err                   --> [1 x 1 x EM_cycles] array with the fractional localization error of the estimate produced at each EM cycle

    % save directory
    save_dir = fullfile('Phase_Error_Analysis','data_out');
    
    % constants
    trials = 50;        % trials per configuration
    mom_samp = 67;      % sample density of aperture plane [Gram-Schmidt] [samples/length]
    pos_samp = 129;     % image plane samples (must be odd!)
    EM_iters_max = 30;  % max EM iterations
    EM_cycles = 15;     % number of times to run the EM algorithm (with different initializations) on the same measurement
    max_order = 5;      % max basis order for GS and Zernike    
    
    % setup apertures
    A = 7;              % area budget   [length]
    D_eff = 20;         % multi-aperture effective diameter [length]
    R_eff = D_eff/2;    % multi-aperture effective radius   [length]
    r = @(n_ap) sqrt(A/pi/n_ap);  % sub-apeture radius                [length]
    
    
    %ap3 = [Polygon(3,0,'radius',R_eff-r(3)),r(3)*ones(3,1)];
    %ring6 = [Polygon(6,0,'radius',R_eff-r(6)),r(6)*ones(6,1)];
    %ring12 = [Polygon(12,0,'radius',R_eff-r(12)),r(12)*ones(12,1)];
    
    mono = [0,0,r(1)];
    golay9 = [Golay9(R_eff-r(9)),r(9)*ones(9,1)];
    ring9 = [Polygon(9,0,'radius',R_eff-r(9)),r(9)*ones(9,1)];
    
    apertures = {mono,golay9,ring9,};
    aperture_names = {'Monolith', 'Golay 9','Ring 9'};     
    
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
    DS.cfg_idx_names = {'Aperture Index (a)','Source Number Index (n)','Min Separation Index (m)','Photon Number Index (p)','Basis Index (b)','Dark Current (e1)','Phase Error (e2)'};
    DS.cfg_data_names = {'Rayleigh Length','Scene','Measurement Mode Counts','Estimated Scene','Log Likelihood','Error'};
    
    
    % array properties (parameter scans)
    DS.num_pho = [1e5,1e6];    % mean photon count
    DS.basis = {'Zernike','Direct-Detection'};
    DS.min_sep_frac = 2.^(linspace(-5,-2,4)); % fractional rayleigh units
    DS.apertures = apertures;                 
    DS.aperture_names = aperture_names;
    DS.num_src = 5;
    DS.dark_lambda = 0;
    DS.phase_sigma = [0,1/100,1/50,1/10]; % [waves]
    DS.cfg_size = [numel(DS.apertures),numel(DS.num_src),numel(DS.min_sep_frac),numel(DS.num_pho),numel(DS.basis),numel(DS.dark_lambda),numel(DS.phase_sigma)]; % the dimensionality of the parameter space range
    DS.data = cell(DS.cfg_size);
    DS.rl = ones(1,numel(DS.apertures))*(2*pi * 1.2197/D_eff); 
    
    
end


