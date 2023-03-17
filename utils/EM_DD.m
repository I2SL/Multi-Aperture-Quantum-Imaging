function [s_xy_trc,s_b_trc,count] = EM_DD(pho_xy, num_sources, aperture, rl, n_em_max, brite_flag, use_gaussian)
% A fast implementation of EM for direct detection.
% The underlying mathematics require that the |psf|^2 has mean equal to mode!
% This allows us to forego discretizing the image plane and finding a
% maximum. We instead just compute the mean of the sources.
% This condition holds true for multi-aperture systems of identical circular sub-apertures
% with radial symmetry and likely holds true for systems of identical
% circular sub-apertures without symmetry.

% initialize source positions
%s_xy = rl/4*(rand(num_sources,2)-.5);
s_xy = 2*rl*(rand(num_sources,2)-.5);

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

end

%{
function [s_b_trc, s_x_trc, s_y_trc, loglike_trc, count] = EM(mode_counts,num_sources,prob_fn,X,Y,rl,n_em_max,brite_flag)
% runs expectation maximization to determine source coordinates and source
% brightnesses from a measurement.
%
% mode_counts       : number of photons detected in each mode
% num_sources       : How many sources involved in the scene
% prob_fn           : modal photon detection PMF given a source located at (x,y)
% X                 : source domain x-coordinate
% Y                 : source domain y-coordinate
% rl                : rayleigh length of the system
% n_em_max          : max number of EM iterations
% brite_flag        : 1 if brightenesses are also to be estimated, 0 otherwise
% ---------------------------------------------
% s_b_trc           : a trace of the brightness estimates across EM iterations
% s_x_trc           : a trace of the source x position estimates across EM iterations
% s_y_trc           : a trace of the source y position estimates across EM iterations
% loglike_trc       : a trace of the log-likelihood of the estimate across EM iterations
% count             : the number of EM iterations performed


% initialize source positions

% radial sub-rayleigh source position initialization
%s_x = rl/4.*cos(2*pi*(0:num_sources-1)/num_sources + pi/2)';
%s_y = rl/4.*sin(2*pi*(0:num_sources-1)/num_sources + pi/2)';


% calculation of log-likelihood offset for multinomial (only necessary for
% logging ln_trc)
sterling = @(n) n.*log(n)-n; % sterling approximation to log(n!)
N = sum(mode_counts);
if N <= 170
    log_N = log(factorial(N));
else
    log_N = sterling(N);
end
log_n = zeros(size(mode_counts));        
log_n(mode_counts <= 170) = log(factorial(mode_counts(mode_counts <= 170)));
log_n(mode_counts > 170) = sterling(mode_counts(mode_counts > 170));
offset = log_N - sum(log_n); % log( N!/ ( n1! ... nM!) )


% random sub-rayleigh source position initialization
s_x = rl/4*(rand(num_sources,1)-.5);
s_y = rl/4*(rand(num_sources,1)-.5);

% initialize source weights
s_b = ones(num_sources,1)/num_sources;

s_x_trc = zeros(num_sources,1);  % source x coordinates
s_y_trc = zeros(num_sources,1);  % source y coordinates
s_b_trc = zeros(num_sources,1);  % source brightnesses
loglike_trc = zeros(1,1);       % log likelihood


% GPU optimization if available
if gpuDeviceCount("available") > 0
    X = gpuArray(X);
    Y = gpuArray(Y);
end



% log probability for all source position
lnP = log(prob_fn(X(:),Y(:)));
lnP = lnP(:,mode_counts >0);

% keep iterating until:
%       - parameters no longer change
%       - the max number of EM iterations have been reached

count = 0;
while ( ~all((s_x - s_x_trc(:,end)) == 0) || ~all((s_y - s_y_trc(:,end)) == 0) ) && count <= n_em_max
    
    % get modal PMFs for each of the source position estimates 
    p_j_est = prob_fn(s_x,s_y); 
    
    % weight the modal PMFS of each source position by the relative source
    % brightness
    p_c = s_b .* p_j_est(:,mode_counts>0); 
    
    % probability per mode given the sources locations
    p_mode = sum(p_c,1);
    
    % log likelihood of the source configuration
    loglike = sum(mode_counts(mode_counts>0).* log(p_mode)) + offset;
    
    % updated traces of the scene parameter estimates accross EM iterations
    s_x_trc(:,count+1) = s_x;                                                  
    s_y_trc(:,count+1) = s_y;                                                 
    s_b_trc(:,count+1) = s_b;                                                  
    loglike_trc(1,count+1) = loglike;  
    
    
    if count < n_em_max    
        % ----------- EXPECTATION STEP -----------

        % measurement weights
        T = mode_counts(mode_counts>0) .* p_c ./ sum(p_c,1);

        % get Q
        Q = sum(pagemtimes(reshape(T,[size(T,1),1,size(T,2)]), reshape(lnP,[1,size(lnP)])),3);

        % reshape Q for 2D    
        Q_2D = zeros([size(X),num_sources]);
        for i = 1:num_sources
            Q_2D(:,:,i) = reshape(Q(i,:),size(X));
        end 

        % ----------- MAXIMIZATION STEP -----------

        %{
        % debugging
        % visualization of objective function for MLE estimators
        peaked_Q = zeros([size(X),num_sources]);
        for i = 1:num_sources
            %ki =  find(Q(i,:) == max(Q(i,:))); % max indices
            peaked_Qi = reshape(Q(i,:),size(X));
            %peaked_Qi(ki) = peaked_Qi(ki) + 1e4;
            peaked_Q(:,:,i) = peaked_Qi; 

            subplot(1,num_sources,i)
            surf(X/rl,Y/rl,peaked_Qi)
            xlabel('x [rl]')
            ylabel('y [rl]')
            zlabel(['Q_',num2str(i)]);
            title(['Intermediate Likelihood for Source ',num2str(i)])
            axis 'square'

        end 
        %}
        
        if brite_flag
            s_b = sum(T,2)/size(T,2); % update source weights
        end
        
        [s_x,s_y] = MLESourceCoords(X,Y,Q_2D);

        nc = size(s_x,3); % number of candidate source locations
        if nc > 1
            % warning('EM iteration recovered degenerate solutions. Selecting first.')
            s_x = s_x(:,1,1);
            s_y = s_y(:,1,1);
        end
        
    end
    
    % increment em updates counter
    count = count + 1;
    
end

if count == n_em_max + 1
    warning('EM reached max number of iterations')
end

s_b_trc = gather(s_b_trc);
s_x_trc = gather(s_x_trc);
s_y_trc = gather(s_y_trc);
loglike_trc = gather(loglike_trc);


end
%}
