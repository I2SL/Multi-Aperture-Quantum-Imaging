clear
addpath('utils/')

% setup apertures
D = 30;     % multi-aperture effective diameter [length]
R = D/2;    % multi-aperture effective radius   [length]
d = 3;      % sub-aperture diameter             [length]
r = 3/2;    % sub-apeture radius                [length]
ap2 = Polygon(2,0,'radius',R);
ap3 = Polygon(3,0,'radius',R);
ap9 = Polygon(9,0,'radius',R);
golay9 = Golay9(R);
aper_coords = {ap2,ap3,ap9,golay9};

% structure for organizing data survey
synth_data = struct();
synth_data.trials = 30;
synth_data.img_samp = 201;
synth_data.subap_samp = 201;
synth_data.subap_radius = r;
synth_data.effap_radius = R;
synth_data.EM_max = 20;
synth_data.max_order = 5;
synth_data.num_pho = [1e3,5e3,1e4,2e4];
synth_data.basis = {'Gram-Schmidt','Zernike','Direct-Detection'};
synth_data.min_sep = 2.^(linspace(-6,-3,8)); % fractional rayleigh units
synth_data.aper_coords = aper_coords;        % [length]
synth_data.num_src = [2,3];
synth_data.data = cell({});



% Run parameter scans using matlab's Parallel Computing Toolbox
visualize = 0;          % visualization trigger
config_index = 1;       % system configuration index (independent of trials)
for b = 1:numel(synth_data.basis)
    for p = 1:numel(synth_data.num_pho)
        for a = 1:numel(synth_data.aper_coords)
            for n = 1:numel(synth_data.num_src)
                for m = 1:numel(synth_data.min_sep)
                    
                    % make the scene
                    n_src = synth_data.num_src(n);
                    src_coords_frac = Polygon(n_src,0,'separation',synth_data.min_sep(m)); % source coordinates in fractional rayleigh units [x/rl]
                    src_brites = ones(n_src,1) / n_src;
                    scene = [src_coords_frac, src_brites];
                    
                    
                    % function inputs
                    np = synth_data.num_pho(p);              % mean photon number                   [integer]
                    mo = synth_data.max_order;               % max modal order                      [integer]
                    bs = synth_data.basis{b};                % basis                                [string] ['Gram-Schmidt','Zernike', 'Direct-Detection']
                    sr = synth_data.subap_radius;            % sub-aperture radius                  [double] [units : length]
                    ac = synth_data.aper_coords{a};          % aperture position                    [Mx2]    [units : length] --> col_1 = kx, col_2 = ky
                    ss = synth_data.subap_samp;              % sub-aperure samples along each axis  --> Bounding box for each aperture has dimensions [subap_samp,subap_samp]
                    is = synth_data.img_samp;                % image-plane samples along each axis  --> Image plane for source position estimates has dimensions [img_samp,img_samp]
                    em = synth_data.EM_max;                  % max number of EM iterations          [integer]
                    vs = visualize;                          % visualization trigger                [boolean]
                    for t = 1:synth_data.trials
                        
                        
                        % run reconstruction 
                        [est_scene, mode_counts, rl, err] = MultiAperture_ConstellationLocalization(...
                            np,...
                            mo,...
                            bs,...
                            sr,...
                            ac,...  
                            scene,...
                            ss,...
                            is,...
                            em,...
                            vs...
                            );
                            
                        
                        % Basis Index, Photon Number Index,  Aperture
                        % Index, Source Number Index, Min Separation Index,
                        % scene, estimated scene, modal measurment, rl,err,
                        % configuration index
                         
                        data(t,:) = {b,p,a,n,m,scene,est_scene,mode_counts,rl,err,config_index};                        
                    end
                    
                    synth_data.data(config_index:config_index-1+synth_data.trials,:) = data;
                    
                    % save current data
                    save('synth_data.mat','synth_data')
                    
                    config_index = config_index + synth_data.trials;
                end
            end
        end
    end
end


% save outputs
    %{
    % save outputs
    save_dir = fullfile(basis,[num2str(ap_num),'-aperture'],[num2str(num_sources),'-src'],[num2str(rl_frac),'rl']);
    mkdir(save_dir)

    meas_file = fullfile(save_dir,'measurement.png');
    est_file = fullfile(save_dir,'estimate.png');
    fig_file = fullfile(save_dir,'figures.fig');
    pst_cnt = 1;
    while isfile(meas_file) || isfile(est_file) || isfile(fig_file)
        meas_file = fullfile(save_dir,['measurement',num2str(pst_cnt),'.png']);
        est_file =  fullfile(save_dir,['estimate',num2str(pst_cnt),'.png']);
        fig_file = fullfile(save_dir,['figures',num2str(pst_cnt),'.fig']);
        pst_cnt = pst_cnt + 1;
    end

    % save figures
    savefig(figs,fig_file)

    % save images
    saveas(figs(1),meas_file);
    saveas(figs(2),est_file);
    %}




% Localization Error v Rayleigh separation (for all modal base and photon
% counts
%load 'synth_data.mat'





% RGYB performance regions - 1 per (aperture configuration, modal basis)
% x - # sources
% y - # min sp
% z - # photons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% survey parameters
%{

pho_mean = 1e3;     % mean photon count for measurement
basis_list = {'Gram-Schmidt','Zernike'};
ap_num_list = [2,3];
src_num_list = [2,3];
rl_frac_list = linspace(1/64,1/16,10);
num_trials = 10;

% Statistical Survey
%load performance_survey_incomplete2.mat
%{
num_entries = num_trials * numel(basis_list) * numel(ap_num_list) * numel(src_num_list) * numel(rl_frac_list);
sz = [num_entries,6];
varTypes = ["string","double","double","double","double","double"];
varNames = ["Basis","# of Apertures","# of Sources","Source Separation (rl)", "# of Photons", "Average Localization Error (rl)"];
T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
%}

k = 1;
for m = 1:numel(basis_list)
    for a = ap_num_list
        for s = src_num_list
            for r = rl_frac_list
                for t = 1:num_trials
                    
                    % static measurement source localization
                    [err,src_coords,est_coords,ap_coords,mode_counts,n_max,rl] = ...
                    MultiAperture_ConstellationLocalization(basis_list{m},a,s,r,pho_mean);
                    
                    % number of photons in measurement
                    pho_count = sum(mode_counts);
                    
                    % add entry to table
                    T(k,:) = {basis_list{m},a,s,r,pho_count,err};
                    k = k+1;
                    
                    
                    save('performance_survey.mat','T');
                    
                end
            end
        end
    end
end 
%}
 


