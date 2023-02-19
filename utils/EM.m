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

end
