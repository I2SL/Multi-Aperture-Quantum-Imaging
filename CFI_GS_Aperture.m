addpath('utils\')

% set interpreter to latex
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


% DARPA specifications
R_max = 10;                 % [m]   maximal radial extent of the multi-aperture configuration
D_max = 2*R_max;            % [m]   maximum diameter of the multi-aperture configuration
A = 7;                      % [m^2] total collection area budget
r = @(n) sqrt(A/n/pi);      % [m]   sub-aperture radius for a given number of apertures(function handle)

% effective rayleigh length
rl = 2*pi * 1.2197/(D_max); % [rads/m] rayleigh length - the rayleigh length is the radial distance (in image space) from the origin to the first zero of the besselj(1) function rl = bsj1_zero/R_eff;   

% constants
max_order = 5;         % Maximum PAD basis order
[nj,mj] = Indices_GramSchmidt(max_order);
num_modes = (max_order+1)^2; % number of modes
mom_samp = 50;    % [samples/m] aperture sampling density
pos_samp = 50;   % [samples / (rad / m)] image plane sampling density


% ------------------------
% Apertures Configurations
% ------------------------
% Each aperture is represented by a [nx3] matrix
% Each row corresponds to a unique sub-aperture
% The columns are organized as [kx-coordinate, ky-coordinate, radius]
mono = [0,0,r(1)];
ring3 = [Polygon(3,0,'radius',R_max-r(3)),r(3)*ones(3,1)];
ring4 = [Polygon(4,0,'radius',R_max-r(4)),r(4)*ones(4,1)];
ring6 = [Polygon(6,0,'radius',R_max-r(6)),r(6)*ones(6,1)];
ring7 = [Polygon(7,0,'radius',R_max-r(7)),r(7)*ones(7,1)];
ring9 = [Polygon(9,0,'radius',R_max-r(9)),r(9)*ones(9,1)];
golay4 = [Golay4(R_max-r(4)),r(4)*ones(4,1)];
golay6 = [Golay6(R_max-r(6)),r(6)*ones(6,1)];
golay7 = [Golay7(R_max-r(7)),r(7)*ones(7,1)];
golay9 = [Golay9(R_max-r(9)),r(9)*ones(9,1)];


apertures = {mono,ring3,ring4,ring6,ring7,ring9,golay4,golay6,golay7,golay9};
aperture_names = {'Monolith','Ring 3','Ring 4','Ring 6','Ring 7','Ring 9','Golay 4','Golay 6','Golay 7','Golay 9'};

% QFI of a circular aperture with radius R_eff
QFI_circ_ap = R_max^2;

% CFI calculation
[X,Y] = meshgrid(rl*linspace(-.5,.5,pos_samp));

%[X,Y] = meshgrid(
%th = linspace(0,2*pi,ngrid);
%r = rl*linspace(.01,.5,ngrid);
%[Th,R] = meshgrid(th,r);
%[X,Y] = polt2cart(Th,R);

%{
CFI = zeros(numel(X),numel(apertures));
for a = 1:numel(apertures)
    [GS_basis_mom,Kx,Ky,d2k,A_tot] = GSbasisMom(max_order,apertures{a},mom_samp);
    s_b = [1;1]/2;
    for i = 1:numel(X)
        alpha_vec = [X(i),Y(i)];
        CFI(i,a) = sum(CFI_r_GramSchmidt(alpha_vec,Kx,Ky,d2k,GS_basis_mom,A_tot,s_b),2);
    end
end


%}

% load in the CFI data
load CFI_GS_aperture.mat

% normalize the CFI by the QFI of the effective aperture size
CFI_QFI = CFI/QFI_circ_ap; % normalized CFI

% make the figure
fig = figure;

% set figure renderer mode
fig.RendererMode = 'manual';
fig.Renderer = 'painters';

% tiled array dimensions
height= 4;
width = 5;
t = tiledlayout(height,width);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% indexing for the tiled array
pos_ap = [1,16,2,3,4,5,17,18,19,20];
pos_cfi= [6,11,7,8,9,10,12,13,14,15];

for a = 1:numel(apertures)    
    
    ax1 = nexttile(pos_ap(a)); 
    VisualizeAperture(apertures{a},R_max);
    set(ax1,'Color','none')
    set(ax1,'YColor','w')
    set(ax1,'XColor','w')

    
    ax2 = nexttile(pos_cfi(a));      
    
    CFI_QFI_a = reshape(CFI_QFI(:,a),size(X));
    axis 'square'
    s = pcolor(X/rl,Y/rl,CFI_QFI_a);
    if 1<= a && a < 3 
        colormap(ax2,winter);
    elseif 3 <= a && a < 7
        colormap(ax2,autumn);
    else
        colormap(ax2,cool);
    end
    s.FaceColor = 'interp';
    s.EdgeColor = 'none';
    xlabel('$x/\sigma_{eff}$')
    ylabel('$y/\sigma_{eff}$')
    axis 'square'
    title(aperture_names{a},'Color','w')
    
    set(ax2,'Color','none')
    set(ax2,'YColor','w')
    set(ax2,'XColor','w')
    
    cbar = colorbar;
    caxis([min(CFI_QFI_a(:)),max(CFI_QFI_a(:))])
    c_rng = [min(CFI_QFI_a(:)),max(CFI_QFI_a(:))];
    c_rng = .9*(c_rng-mean(c_rng)) + mean(c_rng);
    cbar.Ticks = linspace(c_rng(1),c_rng(2),5);
    cbar.TickLabels = arrayfun(@(i) sprintf('%.5f',str2num(cbar.TickLabels{i})),1:numel(cbar.TickLabels),'UniformOutput',false);
    title(cbar,['$\frac{CFI(\Delta)_{',aperture_names{a},'}}{QFI(\Delta)_{Mono_{eff}}}$'],'interpreter','latex','Fontsize',10,'Color','w'); 
    set(cbar,'Color','w')


    
    %c_rng = [min(CFI_QFI_a(:)),max(CFI_QFI_a(:))];
    %c_rng = .9*(c_rng-mean(c_rng)) + mean(c_rng);
    %cbar.Ticks = linspace(c_rng(1),c_rng(2),5);
    %cbar.TickLabels = arrayfun(@(i) sprintf('%.5f',str2num(cbar.TickLabels{i})),1:numel(cbar.TickLabels),'UniformOutput',false);
    
    

    
    
    %{
    
    ax = nexttile;
    CFI_QFI_a = reshape(CFI(:,a)/QFI_circ_ap,size(X));
    s = pcolor(X/rl,Y/rl,CFI_QFI_a);
    s.FaceColor = 'interp';
    s.EdgeColor = 'none';

    
    colormap();
    cbar = colorbar;
    title(cbar,['$\frac{CFI(\Delta)_{',aperture_names{a},'}}{QFI(\Delta)_{mono}}$'],'interpreter','latex','Fontsize',10);
    switch mod(a,5)+1
        case 1
            colormap(ax,spring)
            cbar = colorbar;
            title(cbar,'$\frac{CFI(\Delta)_{mono}}{QFI(\Delta)_{mono}}$','interpreter','latex','Fontsize',10);
        case 2
            colormap(ax,cool)
            cbar = colorbar;
            title(cbar,'$\frac{CFI(\Delta)_{plus 9}}{QFI(\Delta)_{mono}}$','interpreter','latex','Fontsize',10);
        case 3
            colormap(ax,autumn)
            cbar = colorbar;
            title(cbar,'$\frac{CFI(\Delta)_{ring 9}}{QFI(\Delta)_{mono}}$','interpreter','latex','Fontsize',10);
        case 4
            colormap(ax,winter)
            cbar = colorbar;
            title(cbar,'$\frac{CFI(\Delta)_{golay 9}}{QFI(\Delta)_{mono}}$','interpreter','latex','Fontsize',10);
        case 5
            colormap(ax,hsv)
            cbar = colorbar;
            title(cbar,'$\frac{CFI(\Delta)_{ring 3}}{QFI(\Delta)_{mono}}$','interpreter','latex','Fontsize',10);
    end
    c_rng = [min(CFI_QFI_a(:)),max(CFI_QFI_a(:))];
    c_rng = .9*(c_rng-mean(c_rng)) + mean(c_rng);
    cbar.Ticks = linspace(c_rng(1),c_rng(2),5);
    cbar.TickLabels = arrayfun(@(i) sprintf('%.5f',str2num(cbar.TickLabels{i})),1:numel(cbar.TickLabels),'UniformOutput',false);
    %colormap(ax,hsv)


    xlabel('$x/\sigma$')
    ylabel('$y/\sigma$')
    axis 'square'
    %caxis(log([min(CFI(:)),max(CFI(:))]))

    title(aperture_names{a})
    %}
end


set(gcf,'Color','none')
set(gcf, 'InvertHardcopy', 'off');
set(gcf,'Renderer','painter');

title(t,{'PAD Basis CFI for 2-Source Separation'},'Color','w')


function [GS_basis_mom,Kx,Ky,d2k,A_tot] = GSbasisMom(max_order,aperture,mom_samp)
    aper_coords = aperture(:,1:2);
    aper_rads = aperture(:,3);

    % aperture plane discretization
    [Kx,Ky,d2k] = ApertureKxKy(aper_coords,aper_rads,mom_samp);     % Kx [u], Ky [u], d2k [u^2]

    % assign A_tot to be the area of tbe discretized aperture
    A_tot = numel(Kx)*d2k;            

    % Create Gram-Schmidt basis
    GS_basis_mom = genGramSchmidtBasis_mom(max_order,Kx,Ky,d2k);              % basis functions in momentum space
end