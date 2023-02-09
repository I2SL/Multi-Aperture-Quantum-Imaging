function DS = DSformat()
    % DS.cfg_data = {rl, scene, measurement, mle_scene, likelihood, err};
    % cfg = [b,p,a,n,m,t]   --> configuration index vector [b = basis index, p = num_pho index, a = aperture index, n = num_src index, m = min_sep index, t = cfg_trials index]
    % rl                    --> the rayleigh length of the configuration
    % scene                 --> [n_src X 3 ] array contianing the scene parameters : x = scene(:,1), y = scene(:,2), b = scene(:,3) 
    % measurement           --> [1 X n_modes ] array of photon counts in each mode
    % mle_scene             --> [n_src x 3 x EM_cycles] array maximum likelihood estimates arrived at for the EM algorithm 
    % likelihood            --> [1 x 1 x EM_cycles] array with the likelihoods of the estimate produced at each EM cycle
    % err                   --> [1 x 1 x EM_cycles] array with the fractional localization error of the estimate produced at each EM cycle

    
    % constants
    trials = 200;       % trials per configuration
    subap_samp = 101;   % samples per subaperture [Gram-Schmidt] (must be odd!)
    img_samp = 151;     % image plane samples (must be odd!)
    EM_iters = 100;     % max EM iterations
    EM_cycles = 50;     % number of times to run the EM algorithm (with different initializations) on the same 
    max_order = 5;      % max basis order for GS and Zernike       

    % setup apertures
    D_eff = 30;         % multi-aperture effective diameter [length]
    R_eff = D_eff/2;    % multi-aperture effective radius   [length]
    d = 3;      % sub-aperture diameter             [length]
    r = d/2;    % sub-apeture radius                [length]
    ap2 = Polygon(2,0,'radius',R_eff);
    ap3 = Polygon(3,0,'radius',R_eff);
    ap9 = Polygon(9,0,'radius',R_eff);
    golay9 = Golay9(R_eff);
    apertures = {ap3,ap9,golay9};
    aperture_names = {'2 Aperture','3 Aperture','9 Aperture','Golay-9'};


    % data structure with some limited functionality
    DS = struct();
    
    % const properties
    DS.timestamp = datetime;
    DS.trials = trials;
    DS.EM_iters = EM_iters;                 % how many iterations EM should run for
    DS.EM_cycles = EM_cycles;               % how many times we instantiate EM on the same measurement
    DS.img_samp = img_samp;                 % image plane samples (should be odd)
    DS.subap_samp = subap_samp;             % sub aperture samples (should be odd) [for gram-schmidt basis mainly]
    DS.subap_radius = r;                    % sub aperture radius
    DS.ref_unit = r;                        % reference unit
    DS.effap_radius = R_eff;                % effective aperture radius
    DS.max_order = max_order;               % max modal order
    DS.cfg_idx_names = {'Basis Index (b)','Photon Number Index (p)','Aperture Index (a)','Source Number Index (n)','Min Separation Index (m)'};
    DS.cfg_data_names = {'Rayleigh Length','Scene','Measurement Mode Counts','Estimated Scene','Log Likelihood','Error'};
    
    % array properties (parameter scans)
    DS.num_pho = [1e3,5e3,1e4,2e4];
    DS.basis = {'Direct-Detection','Gram-Schmidt','Zernike'};
    DS.min_sep_frac = 2.^(linspace(-6,-3,8)); % fractional rayleigh units
    DS.apertures = apertures;                 % [length]
    DS.aperture_names = aperture_names;
    DS.num_src = 2:5;
    DS.cfg_size = [numel(DS.basis),numel(DS.num_pho),numel(DS.apertures),numel(DS.num_src),numel(DS.min_sep_frac)];  % the dimensionality of the parameter space range
    DS.data = cell(DS.cfg_size);
    
end



%{
function DS_reduced = ReduceDS(DS)
    % Gets the best EM estimate (highest likelihood)
    % for each trial of each configuration
    DS_reduced = DS;

    for i = 1:prod(DS.cfg_size)
        cfg = DS.cfg_data{ind2sub(i,DS.cfg_size)};
        
        for j = 1:size(cfg,1)
            [ml,ml_id] = max(cell2mat(cfg{:,4}),[],3);
            cfg_reduced(j,:) = 
        end
        
        DS_reduced(ind2sub(i,DS.cfg_size)) =
    end
    
end


function [] = measurementAvg(DS)

end


function [] = mergeDS(DS1,DS2)

end


cfg = DS.cfg({b,p,a,n,m}); % returns another cell array with the following properties 
cfg(i,:) = {scene, measurement, mle_scenes, likelihoods, errs};



dd = DS.data{cfg_id{:}};

get estimates with highest loglikelihood among EM cycles
[min_loglike,idx] = min(cell2mat(dd(:,end-1)),[],3)

reduced_data = dd;

for j = 1:DS.EM_cycles
    reduced_data(j,:) = {dd{j,1},dd{j,2},dd{j,3}(:)}
end
DS_reduced = DS;

DS_reduced.data{cfg_id{:}} = data
v = {dd{:,end-1}  }


plot likelihood versus error




% to get the final estimates accross EM cycles we would do do
cfg = DS.cfg({b,p,a,n,m});
[max_like,max_like_id] = max(cell2mat(cfg(i,5)));






DS.cfgs_num = prod(cfg_dim);
DS.cfg_trials = 1:94;            % number of trials for each configuration        
DS.data = cell(DS.cfg_size);

%}







