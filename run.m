%{
basis = 'Zernike'; % basis [Gram-Schmidt, Zernike, Direct-Detection]
ap_num = 2;        % number of apertures 
src_num = 2;       % number of sources
rl_frac = 2/20;   % source spacing (in rayleigh units)
pho_num = 1e3;     % mean photon count for measurement


[src_coords,est_coords,ap_coords,mode_counts,n_max,rl] = ...
MultiAperture_ConstellationLocalization(basis,ap_num,src_num,rl_frac,pho_num);

% compute error
err = LocalizationError(src_coords, est_coords, rl);
%}


%{
num_trials = 10;
basis_list = {'Gram-Schmidt'};
ap_num_list = [2,3];
src_num_list = [2,3];
rl_frac_list = [1/32, 1/64];

for m = 1:numel(basis_list)
    for a = ap_num_list
        for s = src_num_list
            for r = rl_frac_list
                MultiAperture_ConstellationLocalization(basis_list{m},a,s,r,pho_num)
            end
        end
    end
end
%}


% survey parameters
pho_mean = 1e3;     % mean photon count for measurement
%basis_list = {'Gram-Schmidt','Zernike','Direct-Detection'};
basis_list = {'Direct-Detection'};
%ap_num_list = [2,3];
ap_num_list = [3];
src_num_list = [2,3];
rl_frac_list = [1/16, 1/32, 1/64];
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

k = 301;
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

 

