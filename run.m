clear
addpath('utils/')                   % add all of the package functions

D = 30;     % multi-aperture effective diameter [length]
R = D/2;    % multi-aperture effective radius   [length]
d = 3;      % sub-aperture diameter             [length]
r = d/2;    % sub-apeture radius                [length]
aperture = Golay9(R);


% make the scene
n_src = 5;                                  % number of sources
src_brites = ones(n_src,1) / n_src;         % the relative source brightnesses. src_brites should have dimesions [n_src] sum to one and have length n_src.         
src_coords = genMinDistConstellation(src_brites,.2,1); % source coordinates in fractional rayleigh units [x/rl]. The rayleigh length corresponds to the effective diameter
%src_coords = Polygon(5,0,'separation',.1);                     

scene = [src_coords, src_brites];           % input scene. [Nx3] scene(:,1:2)->source coordinates [in fractions of rayleigh], scene(:,3)->relative source brightnesses


% set the system parameters
n_pho = 5e4;        % mean photon number                   [integer]
max_order  = 5;     % max modal order                      [integer]
basis = 'Direct-Detection';  % basis                                [string] ['Gram-Schmidt','Zernike', 'Direct-Detection']
subap_radius = r;   % sub-aperture radius                  [double] [units : length]
aper_coords = aperture; % aperture position                [Mx2]    [units : length] --> col_1 = kx, col_2 = ky
subap_samp = 101;   % sub-aperure samples along each axis  --> Bounding box for each aperture has dimensions [subap_samp,subap_samp]
img_samp = 121;     % image-plane samples along each axis  --> Image plane for source position estimates has dimensions [img_samp,img_samp]
EM_max = 100;        % max number of EM iterations          [integer]
brite_flag = 0;     % brightness estimation flag [boolean]
visualize = 1;      % visualization flag for seeing output [boolean]


[est_scene, mode_counts, rl, err] = MultiAperture_ConstellationLocalization(...
        n_pho,...                   % mean photon number                   [integer]
        max_order,...               % max modal order                      [integer]
        basis,...                   % basis                                [string] ['Gram-Schmidt','Zernike', 'Direct-Detection']
        subap_radius,...            % sub-aperture radius                  [double] [units : length]
        aper_coords,...             % aperture position                    [Mx2]    [units : length] --> col_1 = kx, col_2 = ky   
        scene,...                   % input scene                          [Nx3] scene(:,1:2)->source coordinates [in fractions of rayleigh], scene(:,3)->relative source brightnesses
        subap_samp,...              % sub-aperure samples along each axis  --> Bounding box for each aperture has dimensions [subap_samp,subap_samp]
        img_samp,...                % image-plane samples along each axis  --> Image plane for source position estimates has dimensions [img_samp,img_samp]
        EM_max,...                  % max number of EM iterations          [integer]
        brite_flag,...              % estimate brightness trigger          [boolean]
        visualize...                % visualization trigger                [boolean]
 );