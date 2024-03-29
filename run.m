clear
addpath('utils/')                   % add all of the package functions

% make the aperture
A = 7;      % total aperture area budget        [length^2]
D = 20;     % multi-aperture effective diameter [length]
R = D/2;    % multi-aperture effective radius   [length]
r = @(n_ap) sqrt(A/n_ap/pi); % sub-aperture radius [length]

mono = [0,0,r(1)];
golay9 = [Golay9(R-r(9)),r(9)*ones(9,1)];
ring9 = [Polygon(9,0,'radius',R-r(9)),r(9)*ones(9,1)];
plus9 = [PlusAperture(9,R-r(9)),r(9)*ones(9,1)];
aperture = golay9;

% make the scene
n_src = 2;                                  % number of sources
src_brites = [1e5;1;1;1];                   % the relative source brightnesses. src_brites should have dimesions [n_src] sum to one and have length n_src.         
src_brites = src_brites / sum(src_brites);
delta = 1/4;                               % minimum source separation as a fraction of the rayleigh length
src_coords_frac = genMinDistConstellation(src_brites,delta,1); % source coordinates in fractional rayleigh units [x/rl]. The rayleigh length corresponds to the effective diameter
src_coords_frac = src_coords_frac - src_coords_frac(1,:); 
scene = [src_coords_frac, src_brites];           % input scene. [Nx3] scene(:,1:2)->source coordinates [in fractions of rayleigh], scene(:,3)->relative source brightnesses



% set the system parameters
n_pho = 1e9;        % mean photon number                   [integer]
max_order = 5;      % max modal order                      [integer]
basis = 'Gram-Schmidt';  % basis                           [string] ['Gram-Schmidt','Zernike', 'Direct-Detection']
mom_samp = 67;      % sub-aperure samples along each axis  --> Bounding box for each aperture has dimensions [mom_samp,mom_samp]
pos_samp = 129;     % image-plane samples along each axis  --> Image plane for source position estimates has dimensions [pos_samp,pos_samp]
EM_max = 100;        % max number of EM iterations          [integer]
dark_lambda = 0;    % dark current rate [integer]
phase_sigma = 0;    % phase error (in waves) of mixing channels [double] --> Applies only to local aperture mixing modes 
brite_flag = 0;     % brightness estimation flag [boolean]
exoplanet_flag = 1; % scene is exoplanet search [boolean]
visualize = 1;      % visualization flag for seeing output [boolean]


[est_scene, mode_counts, rl, err] = MultiAperture_ConstellationLocalization(...
        scene,...                   % input scene                          [Nx3] scene(:,1:2)->source coordinates [in fractions of rayleigh], scene(:,3)->relative source brightnesses
        aperture,...                % aperture coordinates and radii       [Mx3] aperture(:,1:2) --> centroid coodrdinates of sub-apertures, aperture(:,3) radius of each sub-aperture    
        n_pho,...                   % mean photon number                   [integer]
        'max_order',max_order,...               % max modal order                      [integer]
        'basis', basis,...                      % basis                                [string] ['Gram-Schmidt','Zernike', 'Direct-Detection']
        'mom_samp', mom_samp,...                % sub-aperure samples along each axis  --> Bounding box for each aperture has dimensions [subap_samp,subap_samp]
        'pos_samp', pos_samp,...                % image-plane samples along each axis  --> Image plane for source position estimates has dimensions [img_samp,img_samp]
        'EM_max', EM_max,...                    % max number of EM iterations          [integer]
        'dark_lambda', dark_lambda,...          % dark current photon arrival rate     [integer]
        'phase_sigma', phase_sigma,...          % phasing error normal random variable's standard deviation [bool]
        'brite_flag', brite_flag,...            % estimate brightness trigger          [boolean]
        'exoplanet_flag', exoplanet_flag,...    % the scene is an exoplanet search     [boolean]
        'visualize', visualize...               % visualization trigger                [boolean]
 );