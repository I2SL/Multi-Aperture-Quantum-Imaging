clear
addpath('utils/')


% constants
trials = 94;         % trials per configuration
subap_samp = 101;   % saples per subaperture
img_samp = 101;     % image plane samples
EM_max = 100;        % max EM iterations
max_order = 5;      % max basis order for GS and Zernike       

% setup apertures

D_eff = 30;     % multi-aperture effective diameter [length]
R_eff = D_eff/2;    % multi-aperture effective radius   [length]
d = 3;      % sub-aperture diameter             [length]
r = d/2;    % sub-apeture radius                [length]
ap2 = Polygon(2,0,'radius',R_eff);
ap3 = Polygon(3,0,'radius',R_eff);
ap9 = Polygon(9,0,'radius',R_eff);
golay9 = Golay9(R_eff);
apertures = {ap2,ap3,ap9,golay9};

% data structure for organizing survey results
DS = struct();
DS.trials = trials;
DS.img_samp = img_samp;
DS.subap_samp = subap_samp;
DS.subap_radius = r;
DS.ref_unit = r;
DS.effap_radius = R_eff;
DS.EM_max = EM_max;
DS.max_order = max_order;
DS.num_pho = [1e3,5e3,1e4,2e4];
DS.basis = {'Gram-Schmidt','Direct-Detection','Zernike'};
DS.min_sep = 2.^(linspace(-6,-3,8)); % fractional rayleigh units
DS.apertures = apertures;        % [length]
DS.num_src = [2,3];
DS.data = cell({});

total_configs = numel(DS.basis)*numel(DS.num_pho)*numel(DS.apertures)*numel(DS.num_src)*numel(DS.min_sep);

parpool(trials)
                   
% (non-dimensionalize) rescale aperture-plane coordinates to the reference unit
subap_radius = DS.subap_radius/DS.ref_unit;   % radius of reference sub-apertures [u]

% Run parameter scans using matlab's Parallel Computing Toolbox
visualize = 0;          % visualization trigger
config_index = 1;       % system configuration index (independent of trials)
for b = 1:numel(DS.basis)
    basis = DS.basis{b};
    for a = 1:numel(DS.apertures)
        
        % make aperture
        ap = DS.apertures{a};
        aper_coords = ap/DS.ref_unit;     % sub-aperture coordinates [u]
        
        % multi-aperture parameters
        ap_num = size(aper_coords,1);           % number of sub-apertures
        A_sub = pi*subap_radius^2;              % subaperture collection area [u^2]
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
        [X,Y] = meshgrid(rl * linspace(-.5,.5,img_samp));       % image plane coordinates [rad/u]

        
        % setup basis
        switch basis
            case 'Gram-Schmidt'

                % indices
                [nj,mj] = Indices_GramSchmidt(max_order);
                                
                % number of modes
                num_modes = numel(nj);
                
                % aperture plane discretization
                [Kx,Ky,d2k] = ApertureKxKy(aper_coords,subap_samp);     % Kx [u], Ky [u], d2k [u^2]
                
                % Create Gram-Schmidt basis
                GS_basis_mom = genGramSchmidtBasis_mom(max_order,Kx,Ky,d2k);                 % basis functions in momentum space
                GS_basis_pos = Basis_GramSchmidt_mom([X(:),Y(:)],Kx,Ky,d2k,GS_basis_mom);    % basis functions in position space
                GS_basis_pos = reshape(GS_basis_pos,[size(X),num_modes]);                    % 2D basis matrix stack

                % probability function
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
        
        
        for n = 1:numel(DS.num_src)
            n_src = DS.num_src(n);
            for m = 1:numel(DS.min_sep)
                    
                    % configure scene
                    src_coords_frac = Polygon(n_src,0,'separation',DS.min_sep(m)); % source coordinates in fractional rayleigh units [x/rl]
                    src_brites = ones(n_src,1) / n_src;
                    scene = [src_coords_frac, src_brites];
                    
                    % source distribution
                    % scene = [s_x,s_y,s_b];
                    src_coords = rl*scene(:,1:2);            % source coordinates [rad/u]
                    num_sources = size(src_coords,1);        % number of sources in the scene
                    s_x = src_coords(:,1); s_y = src_coords(:,2); 
                    s_b = scene(:,3);                        % relative source brightnesses
                                       
                    % make sure min source separation is greater than the resolution of the image plane
                    min_sep = min(pdist(src_coords));
                    dx = X(1,2) - X(1,1);
                    if min_sep < dx
                        warning('Image plane discretization is coarser than minimum source separation')
                    end
                    
                    % get modal probabilities for the given source distribution
                    p_scene = sum(s_b .* prob_fn(s_x,s_y),1);
                    
                    
                for p = 1:numel(DS.num_pho)
                        % configure number of photons
                        n_pho = DS.num_pho(p);
                        
                    disp(['-------Configuration: ' num2str(config_index),'/',num2str(total_configs),'--------'])
                    parfor t=1:trials
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        % simulate the measurement
                        [~, mode_counts] = simulateMeasurement(n_pho, p_scene);

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


                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        % add to data component
                        data(t,:) = {b,p,a,n,m,scene,est_scene,mode_counts,rl,err,config_index};    

                        % display trial completion
                        disp(['Trials Completed: ',num2str(t),'/',num2str(trials)])
                    end
                    
                    % incorporate into data structure
                    DS.data(((config_index-1)*DS.trials+1):(config_index*DS.trials),:) = data;
                    
                    % save current data structure
                    save('Survey_2-3-2023_94MC.mat','DS')
                    
                    % increment the configuration index
                    config_index = config_index + 1;
                end
            end
        end
    end
end