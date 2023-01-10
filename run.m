
modes = 'Direct-Detection'; % modes [Gram-Schmidt, Zernike, Direct-Detection]
ap_num = 1;        % number of apertures 
src_num = 2;       % number of sources
rl_frac = 1/32;    % source spacing (in rayleigh units)
pho_num = 1e4;     % mean photon count for measurement

%MultiAperture_ConstellationLocalization(modes,ap_num,src_num,rl_frac,pho_num) 


mode_list = {'Direct-Detection'};
ap_num_list = [2,3];
src_num_list = [2,3];
rl_frac_list = [1/32, 1/64];

for m = 1:numel(mode_list)
    for a = ap_num_list
        for s = src_num_list
            for r = rl_frac_list
                MultiAperture_ConstellationLocalization(mode_list{m},a,s,r,pho_num)
            end
        end
    end
end
