
addpath('utils/')

% setup apertures
D = 30;     % multi-aperture effective diameter [length]
R = D/2;    % multi-aperture effective radius   [length]
d = 3;      % sub-aperture diameter             [length]
r = 3/2;    % sub-apeture radius                [length]
ap1 = Polygon(1,'radius',R);
ap2 = Polygon(2,0,'radius',R);
ap3 = Polygon(3,0,'radius',R);
ap9 = Polygon(9,0,'radius',R);
golay9 = Golay9(R);
aper_coords = {ap1,ap2,ap3,ap9,golay9};

% structure for organizing data survey
synth_data = struct();
synth_data.trials = 2;
synth_data.img_samp = 201;
synth_data.subap_samp = 201;
synth_data.subap_radius = r;
synth_data.effap_radius = R;
synth_data.EM_max = 20;
synth_data.max_order = 5;
synth_data.num_pho = [1e3,5e3,1e4,2e4];
synth_data.basis = {'Gram-Schmidt','Zernike','Direct-Detection'};
synth_data.min_sep = linspace(1/100, 1/4, 10); % fractional rayleigh units
synth_data.aper_coords = aper_coords;          % [length]
synth_data.num_src = [2,3];
synth_data.data = cell({});


visualize = 1;          % visualization trigger
i = 1;                  % iteration index
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
                    
                    for t = 1:synth_data.trials
                        
                        
                        % run reconstruction 
                        [est_scene, mode_counts, rl, err] = MultiAperture_ConstellationLocalization(...
                            synth_data.num_pho(p),...              % mean photon number                   [integer]
                            synth_data.max_order,...               % max modal order                      [integer]
                            synth_data.basis{b},...                % basis                                [string] ['Gram-Schmidt','Zernike', 'Direct-Detection']
                            synth_data.subap_radius,...            % sub-aperture radius                  [double] [units : length]
                            synth_data.aper_coords{a},...          % aperture position                    [Mx2]    [units : length] --> col_1 = kx, col_2 = ky   
                            scene,...                              % input scene                          [Nx3] scene(:,1:2)->source coordinates [in fractions of rayleigh], scene(:,3)->relative source brightnesses
                            synth_data.subap_samp,...              % sub-aperure samples along each axis  --> Bounding box for each aperture has dimensions [subap_samp,subap_samp]
                            synth_data.img_samp,...                % image-plane samples along each axis  --> Image plane for source position estimates has dimensions [img_samp,img_samp]
                            synth_data.EM_max,...                  % max number of EM iterations          [integer]
                            visualize...                           % visualization trigger                [boolean]
                            );
                        
                        % Basis Index, Photon Number Index,  Aperture
                        % Index, Source Number Index, Min Separation Index,
                        % scene, estimated scene, modal measurment, rl,err,
                        % configuration index
                        synth_data.data(i,:) = {b,p,a,n,m,scene,est_scene,mode_counts,rl,err,config_index};
                        
                        save('synth_data.mat','synth_data')
                    end
                    
                    config_index = config_index + 1;
                end
            end
        end
    end
end




% Localization Error v Rayleigh separation (for all modal base and photon
% counts
load 'synth_data.mat'





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
 


