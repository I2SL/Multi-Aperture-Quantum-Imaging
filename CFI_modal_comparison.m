% Local Zernike Mode vs Gram-Schmidt Modes CFI comparison for 2 point
% source problem using 2 apertures.

clear

% imaging system parameters
subap_radius = 0.5; % radius of sub-apertures [units of length]

% make the sub-aperture length the reference unit
ref_unit = subap_radius;

% rescale everything to the reference unit
subap_radius = subap_radius/ref_unit;   % radius of sub-apertures  [reference unit]
A_sub = pi*subap_radius^2;              % subaperture collection area
rl_1ap =  2*pi*1.22;                                % single-aperture rayleigh length

% multi-aperture parameters
num_apertures = 2;                % number of apertures
A_tot = num_apertures * A_sub;                        % total collection area of the multi-aperture system
U = sqrt(0.5)*[1,1;-1i,1i];                           % aperture mixing matrix

%aper_coords = x(i)*subap_radius*[0,1;0,-1];           % aperture coordinates
%B = pdist(aper_coords);                               % baseline lengths
%max_B = max([1,B]);                                   % max baseline
%rl = rl_1ap/max_B;                                  % multi-aperture rayleigh length



% max modal order
n_max = 5;                                        

% Zernike modes
[n1,m1,v1] = Indices_MixedAperture(n_max, num_apertures);     % mixed mode indices
num_zern_modes = (n_max+1)*(n_max+2)/2 * num_apertures;  % total number of modes


% Gram-Schmidt Modes
[n2,m2] = Indies_GramSchmidt(n_max);
num_GS_modes = (n_max+1)^2;


% Two point source problem parameters
rl_divisor = 10;                                                    % source spacing divisor
alpha_vec = rl_1ap/(2*rl_divisor)*[0,1];                            % source displacement coordinate
n_sources = 2;
s_b = ones(n_sources,1)./ n_sources; % source brightnesses



% loop through aperture separations
N = 50;
x = linspace(1,3,N);
for i = 1:N
% multi-aperture parameters
aper_coords = x(i)*subap_radius*[0,1;0,-1];          % aperture coordinates
B = pdist(aper_coords);                             % baseline lengths
max_B = max([1,B]);                                 % max baseline

% rayleigh length
rl = rl_1ap/max_B;                                  % multi-aperture rayleigh length

% image plane discretization  
ip_dim = 101;
[X,Y] = meshgrid(rl * linspace(-.5,.5,ip_dim));



[Kx,Ky,d2k,GS_basis_mom] = genGramSchmidtBasis2(n_max,aper_coords,301);
GS_basis_pos = Basis_GramSchmidt2(X(:),Y(:),Kx,Ky,d2k,GS_basis_mom);
GS_basis_pos = reshape(GS_basis_pos,[size(X),size(GS_basis_pos,2)]);

%{
% aperture plane discretization
ap_dim = 51;
[aperture,Kx,Ky] = ApertureConfig(aper_coords(:,1),aper_coords(:,2),ap_dim);
rel_ap = subap_radius / Kx(1,end);                      % ratio of sub-aperture radius to aperture plane half-width

[poly_coeff,GS_basis_mom,GS_basis_pos] = genGramSchmidtBasis(Kx,Ky,aperture,n_max,ap_dim,ip_dim,rel_ap);
% ensure the position space modes are properly normalized
% (an on-axis source should produce probability 1 in the 00 mode)
GS_normalization = sqrt(A_tot)*abs(Basis_GramSchmidt(0,0,X,Y,GS_basis_pos(:,:,1)));
GS_basis_pos = GS_basis_pos / GS_normalization;
%}

% Modal contribution to CFI for zernikes
CFI_zern(i,:) = CFI_r_MixedAperture(alpha_vec,n1,m1,v1,U,aper_coords,A_tot,s_b);
CFI_gs(i,:) = CFI_r_GramSchmidt(alpha_vec,X,Y,GS_basis_pos,A_tot,s_b);
end


% Compute QFI
[r,theta] = cart2pol(alpha_vec(:,1),alpha_vec(:,2));
QFI = sum( 4 * ((2*pi)/sqrt(A_sub) *  abs(dr_FTZernike(r,theta,n1,m1))).^2);

% sum CFI over modes
CFI_zern_x = sum(CFI_zern,2);
CFI_gs_x = sum(CFI_gs,2);


% figures
% Show CFI as a function of aperture spacing
figure
plot(x,CFI_zern_x,'red',x,CFI_gs_x,'blue')
title({'Modal Basis CFI Comparison',['Source Separation: $rl_{1ap}/',num2str(rl_divisor),'$']},'interpreter','latex')
xlabel('Aperture Half-Separation [$\delta$ : subaperture radius]','interpreter','latex')
ylabel('Source Half-Separation CFI $\mathcal{J}_{\alpha}$','interpreter','latex')
legend({'Mixed Zernike','Gram-Schmidt'})
ylim([0,(QFI+4*x(end)^2)]);
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

