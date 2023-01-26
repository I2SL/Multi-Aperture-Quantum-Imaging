% Local Zernike Mode vs Gram-Schmidt Modes CFI comparison for 2 point
% source problem using 2 apertures.

% multi-aperture static measurement estimation
nu = 1;             % flux attenuation 

% imaging system parameters
subap_radius = 0.5; % radius of sub-apertures [units of length]

% make the sub-aperture length the reference unit
ref_unit = subap_radius;

% rescale everything to the reference unit
subap_radius = subap_radius/ref_unit;   % radius of sub-apertures  [reference unit]
A_sub = pi*subap_radius^2;              % subaperture collection area

% multi-aperture 
num_apertures = 2;                % number of apertures
A_tot = num_apertures * A_sub;                      % total collection area of the multi-aperture system
U = 0.5*[1,1;-1i,1i];                                 % aperture mixing matrix

% rayleig length
rl_1ap =  2*pi*1.22;                                % single-aperture rayleigh length

% Zernike modes
n_max = 5;                                          % max radial order of zernikes
[n,m,v] = zernikeIndices(n_max, num_apertures);     % mixed mode indices
num_modes = (n_max+1)*(n_max+2)/2 * num_apertures;  % total number of modes

% Two point source problem parameters
%src_coords = rl/20*[0,1;0,-1];
alpha_vec = rl_1ap/20*[0,1];                            % source displacement coordinate
n_sources = 2;
b_s = ones(n_sources,1)./ n_sources; % source brightnesses

% loop through aperture separations
N = 100;
x = linspace(.001,1,N);
for i = 1:N
% multi-aperture parameters
aper_coords = x(i)*subap_radius*[0,1;0,-1];            % aperture coordinates
B = pdist(aper_coords);                             % baseline lengths
max_B = max([1,B]);                                 % max baseline

% rayleigh length
rl = rl_1ap/max_B;                                  % multi-aperture rayleigh length

% Modal contribution to CFI for zernikes
CFI_zern(i,:) = MixedApertureCFI(alpha_vec,n,m,v,U,aper_coords,A_tot,b_s);
end

CFI_x = sum(CFI_zern,2);

% Compute QFI
[n1,m1,~] = zernikeIndices(n_max, 1);     % single aperture indices
[r,theta] = cart2pol(alpha_vec(:,1),alpha_vec(:,2));
QFI = sum( 4 * ((2*pi)/sqrt(A_sub) *  abs(dr_FTZernikeBasis(r,theta,n1,m1))).^2);

k1 = 1;
k2 = num_apertures;
for i=1:num_modes/num_apertures
CFI_zern_nm(:,i) = sum(CFI_zern(:,k1:k2),2);
k1 = k2+1;
k2 = k1 + num_apertures - 1;
end


% figures


% Show CFI as a function of aperture spacing
figure
plot(x,CFI_x)
title('2-Aperture Local-Mode CFI')
xlabel('Aperture Half Spacing [$\delta$ : subaperture radius]','interpreter','latex')
ylabel('Source Half Separation CFI $\mathcal{J}_{\alpha \alpha}$','interpreter','latex')
ylim([0,0.5*(QFI+4*x(end))]);


% show CFI per mode at different aperture separations
figure
bar3(x,CFI_zern_nm(:,1:6))
axis 'square'
zlabel('Source Half Separation CFI by Mode $[\mathcal{J}_{\alpha \alpha}]_j$','interpreter','latex')
ylabel('Aperture Half Spacing [$\delta$ : subaperture radius]','interpreter','latex')
xlabel('Local Mode Index $j$','interpreter','latex')
xticklabels(0:5)
title('Aperture-Sum CFI by Local Mode')







function CFI_r_nm = GramSchmidtCFI(x,y,X,Y,GS_modes,A_tot,b_s)
    % GS correlation function  (we have this)
    p_nm = sum(b_s .* GramSchmidtModalProb(x,y,X,Y,GS_modes,A_tot));
    
    % Take the numerical derivative of the probability wrt to the radial coordinate
    [theta,~] = cart2pol(x,y);
    dr = 1e-2; % some small parameter
    dx = dr*cos(theta);
    dy = dr*sin(theta);
    dr_p_nm_s = (GramSchmidtModalProb(x_dx,y_dy,X,Y,GS_modes,A_tot) - p_nm) / dr;
    dr_p_nm = sum(b_s.*dr_p_nm_s);
    
    CFI_r_nm = (dr_p_nm).^2 ./ p_nm;
    
    
end







function p = GramSchmidtModalProb(x,y,X,Y,GS_modes,A_tot)
    % Computes the modal probability of detecting a photon in the
    % Gram-Schmidt basis for sources positioned at [x,y].
    %
    % x,y - source locations
    % X,Y - image plane meshgrids
    % GS_modes  - a matrix stack representing the GS modes over X,Y
    % A_tot - total area of the aperture
    
    dx = X(1,2) - X(1,1);
    dy = Y(2,1) - Y(1,1);

    n_modes = size(GS_modes,3);
        
    p = zeros([numel(x),n_modes]);
    for i = 1:n_modes
        GS_i = GS_modes(:,:,i);
        p(:,i) = interp2(X,Y, (2*pi / sqrt(A_tot) * abs(GS_i)).^2,x,y) * dx * dy;
    end
    
    p = p ./ sum(p,2);              % normalize
end



function CFI_r_nm_v = MixedApertureCFI(alpha_vec,n,m,v,U,aper_coords,A_tot,b_s)

    
    % Probability of mixed mode (n,m,v)
    P_k_nu = sum(b_s .* abs(corrFnMixedApertureBasis(alpha_vec,n,m,v,U,aper_coords,A_tot)).^2);
    
    % partial derivative of the probability with respect to the separation
    
    dr_P_k_nu_s = real(   conj(corrFnMixedApertureBasis(alpha_vec,n,m,v,U,aper_coords,A_tot)) .* dr_corrFnMAB(alpha_vec,n,m,v,U,aper_coords,A_tot) +...
                        conj(corrFnMixedApertureBasis(-alpha_vec,n,m,v,U,aper_coords,A_tot)).* dr_corrFnMAB(-alpha_vec,n,m,v,U,aper_coords,A_tot)...
                     );
            
    % weight by the source brightnesses
    %dr_P_k_nu = sum(b_s .* dr_P_k_nu_s);
    dr_P_k_nu = dr_P_k_nu_s; 
    
    % 2-point source separation CFI
    CFI_r_nm_v = (dr_P_k_nu).^2 ./ P_k_nu;

end


function Gamma_nm_v = corrFnMixedApertureBasis(xy_coords,n,m,v,U,aper_coords,A_tot)
       Gamma_nm_v =  2*pi/sqrt(A_tot) * MixedApertureBasis(xy_coords,n,m,v,U,aper_coords);
end


function dr_Gamma_nm_v = dr_corrFnMAB(xy_coords,n,m,v,U,aper_coords,A_tot)
% derivative of Gamma wrt to source half-separation
[theta,r] = cart2pol(xy_coords(:,1),xy_coords(:,2));
dr_Gamma_nm_v =  2*pi/sqrt(A_tot) * conj(dr_FTZernikeBasis(r,theta,n,m) .* phaseFn(xy_coords,v,U,aper_coords)...
                                     + FTZernikeBasis(r,theta,n,m) .* dr_phaseFn(xy_coords,v,U,aper_coords));    
end


function psi_nm_v = MixedApertureBasis(xy_coords,n,m,v,U,aper_coords)
    % xy_coords : Nx2 matrix with cartesian coordinate pairs at which to
    % evaluate the basis functions
    % n,m : 1xM vectors containing Zernike mode indices
    % v   : 1xM vectors containing mixed mode indices 
    % U   : KxK unitary mixing matrix defining how the local modes will be mixed
    % aper_coords: Kx2 matrix with cartesian coordinate pairs for the aperture positions
    
    % phase from multi-aperture mixing
    B = phaseFn(xy_coords,v,U,aper_coords);
    
    % mixed mode
    [theta,r] = cart2pol(xy_coords(:,1),xy_coords(:,2));
    psi_nm_v = B .* FTZernikeBasis(r,theta,n,m);
end

function B = phaseFn(xy_coords,v,U,aper_coords)

    % phase introduced by each shifted sub-aperture
    phase = exp( 1i * aper_coords * xy_coords');
    
    % aperture mixing
    B = (U(v,:) * phase).';
end

function dr_B = dr_phaseFn(xy_coords,v,U,aper_coords)    
    % phase introduced by each shifted sub-aperture
    phase = exp( 1i * aper_coords * xy_coords');
    
    % direction projections
    xy_hat_coords = xy_coords./vecnorm(xy_coords,2,2);
    aper_cos = 1i * aper_coords * xy_hat_coords';
 
    % derivative of modal interference term with respect to source separation
    dr_B = (U(v,:) * (aper_cos .* phase)).';
end

function dr_z = dr_FTZernikeBasis(r,theta,n,m)
    % derivative of the zernike basis modes with respect to the radial
    % parameter
    dr_z = ((-1).^(n/2)) ./(4*sqrt(pi)*sqrt(n+1)) .* FTzAngle(theta,m) .* (besselj(n-1,r) - besselj(n+3,r));
end


function z = FTZernikeBasis(r,theta,n,m)
    % returns the evaluation of the function and its linear indices
    %
    % INPUTS
    % n - radial index          (row vector 1xK)
    % m - azimuthal index       (row vector 1xK)
    % r - radial argument       (col vector Dx1)
    % theta - angular argument  (col vector Dx1)
    % (r,theta) are coordinate pairs, (n,m) are index pairs
    %
    % OUTPUTS
    % z - the function output   (matrix DxK)
    
    %z = radius*FTzRadial(radius*r,n) .* FTzAngle(theta,m);    
    z = FTzRadial(r,n) .* FTzAngle(theta,m);    
end

function u = FTzRadial(r,n)
    
    % sinc-bessel in polar
    J = besselj(repmat(n+1,[size(r,1),1]), repmat(r,[1,size(n,2)]))./ repmat(r,[1,size(n,2)]);
    
    % fill in singularities
    J(r==0,n+1 == 1) = 0.5;
    J(r==0,n+1 > 1) = 0;
    
    % radial function
    u = (-1).^(n./2) .* sqrt((n+1)/pi) .*  J;
end

function v = FTzAngle(theta,m)
    v = zeros(numel(theta),numel(m));
    
    % angular function
    c = cos(abs(m).*theta);
    s = sin(abs(m).*theta);
    
    v(:,m>0) = sqrt(2) * c(:,m>0);
    v(:,m==0)= 1;
    v(:,m<0) = sqrt(2) * s(:,m<0);
end

function [nj,mj,vj] = zernikeIndices(n_max,n_apertures)
    nj = [];
    mj = [];

    for n = 0:n_max
        for m = -n:2:n
            nj = [nj, n];
            mj = [mj, m];
        end
    end

    n_modes = numel(nj);

    % radial, azimuthal, and aperture index lists
    [vj,jj] = ndgrid(1:n_apertures,1:n_modes);
    nj = nj(jj(:));
    mj = mj(jj(:));
    vj = vj(:)';
    
end
