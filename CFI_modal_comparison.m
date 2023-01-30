% Local Zernike Mode vs Gram-Schmidt Modes CFI comparison for 2 point
% source problem using 2 apertures.

clear
addpath('utils/')

% imaging system parameters
subap_radius = 0.5; % radius of sub-apertures [units of length]

% make the sub-aperture length the reference unit
ref_unit = subap_radius; % the reference unit [u]

% rescale everything to the reference unit
subap_radius = subap_radius/ref_unit;   % radius of sub-apertures [u]
A_sub = pi*subap_radius^2;              % subaperture collection area [u^2]
rl_1ap =  2*pi*1.22;                    % single-aperture rayleigh length [rad/u]

% multi-aperture parameters
num_apertures = 2;                      % number of apertures
A_tot = num_apertures * A_sub;          % total collection area of the multi-aperture system
U = sqrt(0.5)*[1,1;-1i,1i];             % aperture mixing matrix

% max modal order
n_max = 5;                                        

% Zernike modes
[n1,m1,v1] = Indices_MixedAperture(n_max, num_apertures);       % mode indices
num_zern_modes = (n_max+1)*(n_max+2)/2 * num_apertures;         % total number of modes

% Gram-Schmidt Modes
[n2,m2] = Indices_GramSchmidt(n_max);                            % mode indices
num_GS_modes = (n_max+1)^2;                                     % total number of modes

% Two point source problem parameters
n_sources = 2;
rl_divisor = 10;                                                    % source spacing divisor
alpha_vec = rl_1ap/(2*rl_divisor)*[0,1];                            % source displacement coordinate
s_b = ones(n_sources,1)./ n_sources; % source brightnesses

% loop through aperture separations
num_spacings = 50;
x = linspace(1,3,num_spacings); % list of aperture half-separations
CFI_zern = zeros(num_spacings, num_zern_modes);
CFI_gs = zeros(num_spacings, num_GS_modes);
for i = 1:num_spacings
    
% multi-aperture parameters
aper_coords = x(i)*subap_radius*[0,1;0,-1];         % aperture coordinates
B = pdist(aper_coords);                             % baseline lengths
max_B = max([1,B]);                                 % max baseline

% aperture plane discretization
subap_sampling = 301;     % number of k-space samples across each circular hard sub-aperture (total samples is subap_sampling^2)
[kx,ky,d2k] = ApertureKxKy(aper_coords,subap_sampling); 


% Modal contribution to CFI for mixed Zernike modes
CFI_zern(i,:) = CFI_r_MixedAperture(alpha_vec,n1,m1,v1,U,aper_coords,A_tot,s_b);

% Modal contribution to CFI for Gram-Schmidt modes
GS_basis_mom = genGramSchmidtBasis_mom(n_max,kx,ky,d2k); % generate GS basis
CFI_gs(i,:) = CFI_r_GramSchmidt(alpha_vec,kx,ky,d2k,GS_basis_mom,A_tot,s_b);
end


% Compute single-aperture QFI
[r,theta] = cart2pol(alpha_vec(:,1),alpha_vec(:,2));
QFI_1ap = sum( 4 * ((2*pi)/sqrt(A_sub) *  abs(dr_FTZernike(r,theta,n1,m1))).^2);

% sum CFI over modes
CFI_zern_x = sum(CFI_zern,2);
CFI_gs_x = sum(CFI_gs,2);


% figures
% Show CFI as a function of aperture spacing
figure
plot(x,CFI_gs_x,'b',x,CFI_zern_x,'--r','LineWidth',2)
title({'Modal Basis CFI Comparison',['Source Separation: $rl_{1ap}/',num2str(rl_divisor),'$']},'interpreter','latex')
xlabel('Aperture Half-Separation [$\delta$ : subaperture radius]','interpreter','latex')
ylabel('Source Half-Separation CFI $\mathcal{J}_{\alpha}$','interpreter','latex')
legend({'Gram-Schmidt','Mixed Zernike'})
ylim([0,(QFI_1ap+4*x(end)^2)]);
xlim([min(x),max(x)])


% show CFI per mode at different aperture separations
figure
subplot(1,2,1)
[~,ib]=ismember([0,0,1;0,0,2;1,-1,1;1,-1,2],[n1;m1;v1].','rows');
bar3(x,CFI_zern(:,ib))
axis 'square'
title({'CFI by Mode','Mixed Zernike Local Modes'},'interpreter','latex')
xlabel('Mixed Local Mode Index $(n_z,m_z,\nu)$','interpreter','latex')
ylabel('Aperture Half-Separation [$\delta$ : subaperture radius]','interpreter','latex')
zlabel('Source Half-Separation CFI by Mode $[\mathcal{J}_{\alpha}]_{n_z,m_z,\nu}$','interpreter','latex')
xticklabels({'(0,0),+','(0,0),-','(0,-1),+','(0,-1),-'})
ylim([min(x),max(x)])

subplot(1,2,2)
[~,ib2]=ismember([0,0;0,1;0,2;0,3],[n2;m2].','rows');
bar3(x,CFI_gs(:,ib2))
axis 'square'
title({'CFI by Mode','Gram-Schmidt PAD Modes'},'interpreter','latex')
xlabel('Mixed Local Mode Index $(n,m)$','interpreter','latex')
ylabel('Aperture Half-Separation [$\delta$ : subaperture radius]','interpreter','latex')
zlabel('Source Half-Separation CFI by Mode $[\mathcal{J}_{\alpha}]_{n,m}$','interpreter','latex')
xticklabels({'(0,0)','(0,1)','(0,2)','(0,3)'})
ylim([min(x),max(x)])

