clear
addpath('utils/')                   % add all of the package functions


% make the aperture
D = 30;     % multi-aperture effective diameter [length]
R = D/2;    % multi-aperture effective radius   [length]
d = 3;      % sub-aperture diameter             [length]
r = d/2;    % sub-apeture radius                [length]
golay9 = [Golay9(R-r),r*ones(9,1)];    
aperture = golay9;

% make the scene
n_src = 3;                                  % number of sources
src_brites = ones(n_src,1) / n_src;         % the relative source brightnesses. src_brites should have dimesions [n_src] sum to one and have length n_src.         
src_coords_frac = genMinDistConstellation(src_brites,1/16,1); % source coordinates in fractional rayleigh units [x/rl]. The rayleigh length corresponds to the effective diameter
scene = [src_coords_frac, src_brites];           % input scene. [Nx3] scene(:,1:2)->source coordinates [in fractions of rayleigh], scene(:,3)->relative source brightnesses


% set the system parameters
n_pho = 2e3;        % mean photon number                   [integer]
max_order  = 5;     % max modal order                      [integer]
basis = 'Direct-Detection';  % basis                           [string] ['Gram-Schmidt','Zernike', 'Direct-Detection']
mom_samp = 67;   % sub-aperure samples along each axis  --> Bounding box for each aperture has dimensions [subap_samp,subap_samp]
pos_samp = 129;     % image-plane samples along each axis  --> Image plane for source position estimates has dimensions [img_samp,img_samp]
EM_max = 100;        % max number of EM iterations          [integer]
brite_flag = 0;     % brightness estimation flag [boolean]
visualize = 1;      % visualization flag for seeing output [boolean]


[est_scene, mode_counts, rl, err] = MultiAperture_ConstellationLocalization(...
        n_pho,...                   % mean photon number                   [integer]
        max_order,...               % max modal order                      [integer]
        basis,...                   % basis                                [string] ['Gram-Schmidt','Zernike', 'Direct-Detection']
        aperture,...                % aperture coordinates and radii       [Mx3] aperture(:,1:2) --> centroid coodrdinates of sub-apertures, aperture(:,3) radius of each sub-aperture    
        scene,...                   % input scene                          [Nx3] scene(:,1:2)->source coordinates [in fractions of rayleigh], scene(:,3)->relative source brightnesses
        mom_samp,...              % sub-aperure samples along each axis  --> Bounding box for each aperture has dimensions [subap_samp,subap_samp]
        pos_samp,...                % image-plane samples along each axis  --> Image plane for source position estimates has dimensions [img_samp,img_samp]
        EM_max,...                  % max number of EM iterations          [integer]
        brite_flag,...              % estimate brightness trigger          [boolean]
        visualize...                % visualization trigger                [boolean]
 );