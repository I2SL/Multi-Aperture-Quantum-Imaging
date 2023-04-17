function [s_b_trc, s_x_trc, s_y_trc, loglike_trc, count] = EM_DD(pho_xy, num_sources, aperture, rl, n_em_max, brite_flag)
% A fast implementation of EM for direct detection.
% The underlying mathematics require that the |psf|^2 has mean equal to
% mode. This allows us to forego discretizing the image plane and finding a
% maximum. We instead just compute the mean of the photon arrival locations.
% 
% (This condition holds true for multi-aperture systems of identical circular sub-apertures
% with radial symmetry and likely holds true for systems of identical
% circular sub-apertures without symmetry.)

% ------------------------
% INPUTS:
% ------------------------
% pho_xy        : Nx2 matrix of photon arrival locations at the image plane
% num_sources   : number of sources in the scene
% aperture      : the system aperture configuration
% n_em_max      : maximum number of EM iterations allowed
% brite_flag    : a boolean indicating whether brightnesses are to be estimated


% initialize source positions
s_xy = rl/4*(rand(num_sources,2)-.5);

% initialize source weights
s_b = ones(num_sources,1)/num_sources;

% keep track of source coordinates, brightnesses, and log likelihood
s_xy_trc = zeros(num_sources,2);    % source xy coordinates
s_b_trc = zeros(num_sources,1);     % source brightnesses
loglike_trc = zeros(1);             % loglikelihood of constellation estimate

% GPU optimization if available
if gpuDeviceCount("available") > 0
    pho_xy = gpuArray(pho_xy);
end

count = 0;
while (~all(s_xy - s_xy_trc(:,:,end)==0,'all')) && count <= n_em_max
    
    % track of estimate
    s_b_trc(:,count+1) = s_b;
    s_xy_trc(:,:,count+1) = s_xy;

    
    % E-step
    P_rj = zeros(num_sources, size(pho_xy,1));              % probability of photon arrival position r given it was emitted from source j
    for j = 1:num_sources
        P_rj(j,:) = abs(MultiAperturePSF(pho_xy - s_xy(j,:),aperture).').^2;
        
    end
    weights = (P_rj .* s_b) ./ sum(P_rj .* s_b , 1);        % source assignment probabilities. Each column of the weights matrix is the probability distribution over the source origin for a photon
    
    
    % save log likelihood of configuration     
    loglike_trc(count+1) = sum(log(sum(P_rj .* s_b,1)),2); % + const.
    
    % M-step
    if brite_flag
        s_b = mean(weights,2);
    end
    s_xy = squeeze(sum(weights .* permute(pho_xy,[3,1,2]), 2) ./ sum(weights,2));

    % increment iteration count
    count = count + 1;   
        
end

if count == n_em_max + 1
    warning('EM reached max number of iterations')
end

% retreive outputs from gpu 
s_b_trc = gather(s_b_trc);
s_xy_trc = gather(s_xy_trc);
s_x_trc = squeeze(s_xy_trc(:,1,:));
s_y_trc = squeeze(s_xy_trc(:,2,:));
end


% Same as above but creates a live video of the estimation process
%{
function [s_b_trc, s_x_trc, s_y_trc, count] = EM_DD(pho_xy, num_sources, aperture, rl, n_em_max, brite_flag, use_gaussian)
% A fast implementation of EM for direct detection.
% The underlying mathematics require that the |psf|^2 has mean equal to mode!
% This allows us to forego discretizing the image plane and finding a
% maximum. We instead just compute the mean of the photon arrival locations.
% This condition holds true for multi-aperture systems of identical circular sub-apertures
% with radial symmetry and likely holds true for systems of identical
% circular sub-apertures without symmetry.



% initialize source positions
s_xy = rl/4*(rand(num_sources,2)-.5);

% initialize source weights
s_b = ones(num_sources,1)/num_sources;

% keep track of source coordinate, brightnesses
s_xy_trc = zeros(num_sources,2);    % source xy coordinates
s_b_trc = zeros(num_sources,1);     % source brightnesses

% GPU optimization if available
if gpuDeviceCount("available") > 0
    pho_xy = gpuArray(pho_xy);
end

if use_gaussian
    % gaussian psf
    psf = @(xy) sqrt(mvnpdf(xy,[],rl^2*eye(2)));
else
    % point spread function
    psf = @(xy) MultiAperturePSF(xy,aperture);
end

%{
v = VideoWriter('EM_2src_DD.avi');
v.FrameRate = 2;
open(v)
%}

count = 0;
while (~all(s_xy - s_xy_trc(:,:,end) == 0,'all')) && count <= n_em_max
    
    
    % keep track of estimates per iteration
    s_b_trc(:,count+1) = s_b;
    s_xy_trc(:,:,count+1) = s_xy;
    

    % E-step
    P_rj = zeros(num_sources, size(pho_xy,1));              % probability of photon arrival position r given it was emitted from source j
    for j = 1:num_sources
        P_rj(j,:) = abs(psf(pho_xy - s_xy(j,:)).').^2;
    end
    weights = (P_rj .* s_b) ./ sum(P_rj .* s_b , 1);        % source assignment probabilities. Each column of the weights matrix is the probability distribution over the source origin for a photon

    %{
    % visualize source weights (for 2 source only)
    figure;
    red = [1,0,0];
    blue = [0,0,1];
    colormap([linspace(red(1),blue(1),1000)',linspace(red(2),blue(2),1000)',linspace(red(3),blue(3),1000)']);
    cbar = colorbar;
    cbar.XTick = 0:1;

    hold on
    scatter(pho_xy(:,1)/rl,pho_xy(:,2)/rl,1,weights(1,:)'.*red+weights(2,:)'.*blue)
    scatter(s_xy(:,1)/rl,s_xy(:,2)/rl,'filled','black')
    hold off
    xlabel('x [rl]')
    ylabel('y [rl]')
    title({'EM on Direct Detection',['Sources: 2'],['Iteration: ' num2str(count)]})
    axis square
    ylabel(cbar,['$P^{[',num2str(count),']}(s=1|\vec{r}^{(i)})$'],'interpreter','latex')
        
    frame = getframe(gcf);
    writeVideo(v,frame)
    close all
    %}
    
    
    % M-step
    if brite_flag
        s_b = mean(weights,2);
    end
    s_xy = squeeze(sum(weights .* permute(pho_xy,[3,1,2]), 2) ./ sum(weights,2));

    % increment iteration count
    count = count + 1;   


end

%{
close(v)
%}

if count == n_em_max + 1
    warning('EM reached max number of iterations')
end

s_b_trc = gather(s_b_trc);
s_xy_trc = gather(s_xy_trc);
s_x_trc = squeeze(s_xy_trc(:,1,:));
s_y_trc = squeeze(s_xy_trc(:,2,:));
end
%}




