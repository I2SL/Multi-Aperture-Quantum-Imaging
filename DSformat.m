function DS = DSformat()
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
    save_dir = fullfile('Survey_4-16-23_9ap_configs_src_exitance_Mono_addendum','data_out');
    
    % constants
    trials = 50;        % trials per configuration
    mom_samp = 67;      % sample density of aperture plane [Gram-Schmidt] [samples/length]
    pos_samp = 129;     % image plane samples (must be odd!)
    EM_iters_max = 30;  % max EM iterations
    EM_cycles = 15;     % number of times to run the EM algorithm (with different initializations) on the same measurement
    max_order = 5;      % max basis order for GS and Zernike    
    
    % setup apertures
    D_eff = 30;         % multi-aperture effective diameter [length]
    R_eff = D_eff/2;    % multi-aperture effective radius   [length]
    d = 3;              % sub-aperture diameter             [length]
    r = d/2;            % sub-apeture radius                [length]
    
    
    ap3 = [Polygon(3,0,'radius',R_eff-r),r*ones(3,1)];
    ring6 = [Polygon(6,0,'radius',R_eff-r),r*ones(6,1)];
    ring12 = [Polygon(12,0,'radius',R_eff-r),r*ones(12,1)];
    
    mono = [0,0,r];
    plus9 = [PlusAperture(9,R_eff-r),r*ones(9,1)];
    ring9 = [Polygon(9,0,'radius',R_eff-r),r*ones(9,1)];
    golay9 = [Golay9(R_eff-r),r*ones(9,1)];
    
    
    apertures = {mono};
    aperture_names = {'mono_9PhoFlux'};    
    
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
    DS.cfg_idx_names = {'Aperture Index (a)','Source Number Index (n)','Min Separation Index (m)','Photon Number Index (p)','Basis Index (b)'};
    DS.cfg_data_names = {'Rayleigh Length','Scene','Measurement Mode Counts','Estimated Scene','Log Likelihood','Error'};
    
    
    % array properties (parameter scans)
    DS.num_pho = 0.0900*[1e4,1e5,1e6];    % mean photon count
    DS.basis = {'Gram-Schmidt','Direct-Detection'};
    DS.min_sep_frac = 2.^(linspace(-6,-3,7)); % fractional rayleigh units
    DS.apertures = apertures;                 
    DS.aperture_names = aperture_names;
    DS.num_src = 3:5;
    DS.cfg_size = [numel(DS.apertures),numel(DS.num_src),numel(DS.min_sep_frac),numel(DS.num_pho),numel(DS.basis)]; % the dimensionality of the parameter space range
    DS.data = cell(DS.cfg_size);
    
   % construct rayleigh lengths for all apertures and store in a field
   DS.rl = zeros(1,numel(DS.apertures)); % rayleigh lengths
   DS.D_eff = zeros(1,numel(DS.apertures)); % effective aperture diameters
   
   % get the effective aperture diameter
   for a = 1:numel(DS.apertures)
        ap = DS.apertures{a};
        ap_num = size(ap,1);
        aper_coords = ap(:,1:2);
        aper_rads = ap(:,3);
        
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
        DS.rl(a) = 2*pi * 1.2197/D_eff; % rayleigh length in units of [rads/length] 
        DS.D_eff(a) = D_eff;
   end
   
    
    
    
    
    
end


