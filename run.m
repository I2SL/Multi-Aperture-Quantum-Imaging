clear
addpath('utils/')                   % add all of the package functions

% make the aperture
D = 30;     % multi-aperture effective diameter [length]
R = D/2;    % multi-aperture effective radius   [length]
d = 3;      % sub-aperture diameter             [length]
r = d/2;    % sub-apeture radius                [length]
golay9 = [Golay9(R-r),r*ones(9,1)];
ap9 = [Polygon(9,0,'radius',R-r),r*ones(9,1)];
ap3 = [Polygon(3,0,'radius',R-r),r*ones(3,1)];
ap7 = [Polygon(7,0,'radius',R-r),r*ones(7,1)];
ap5 = [Polygon(5,0,'radius',R-r),r*ones(5,1)];
plus9 = [PlusAperture(9,R-r),r*ones(9,1)];
aperture = plus9;

% make the scene

n_src = 5;                                  % number of sources
src_brites = ones(n_src,1) / n_src;         % the relative source brightnesses. src_brites should have dimesions [n_src] sum to one and have length n_src.         
delta = 1/32;                               % minimum source separation as a fraction of the rayleigh length
src_coords_frac = genMinDistConstellation(src_brites,delta,1); % source coordinates in fractional rayleigh units [x/rl]. The rayleigh length corresponds to the effective diameter
scene = [src_coords_frac, src_brites];           % input scene. [Nx3] scene(:,1:2)->source coordinates [in fractions of rayleigh], scene(:,3)->relative source brightnesses

%load('test_scene.mat');

% set the system parameters
n_pho = 6e5;        % mean photon number                   [integer]
max_order = 1;      % max modal order                      [integer]
basis = 'Gram-Schmidt';  % basis                           [string] ['Gram-Schmidt','Zernike', 'Direct-Detection']
%basis = 'Direct-Detection';  % basis                           [string] ['Gram-Schmidt','Zernike', 'Direct-Detection']
mom_samp = 67;   % sub-aperure samples along each axis  --> Bounding box for each aperture has dimensions [mom_samp,mom_samp]
pos_samp = 129;     % image-plane samples along each axis  --> Image plane for source position estimates has dimensions [pos_samp,pos_samp]
EM_max = 30;        % max number of EM iterations          [integer]
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

save('est_scene_DD_golay9.mat','est_scene');