function [s_b_trc, s_x_trc, s_y_trc, count] = EM(mode_counts,num_sources,prob_fn,X,Y,rl,n_em_max)
% runs expectation maximization to determine source coordinates and source
% brightnesses from a measurement.
%
% mode_counts       : 
% num_sources       : How many sources involved in the scene
% prob_fn           : modal photon detection PMF given a source located at (x,y)
% X                 :
% Y                 :
% rl                :
% n_em_max          : max number of EM iterations

% initialize source positions

% random sub-rayleigh source position initialization
s_x = rl*(rand(num_sources,1)-.5);
s_y = rl*(rand(num_sources,1)-.5);

% radial sub-rayleigh source position initialization
%s_x = rl/4.*cos(2*pi*(0:num_sources-1)/num_sources)';
%s_y = rl/4.*sin(2*pi*(0:num_sources-1)/num_sources)';

% initialize source weights
s_b = ones(num_sources,1)/num_sources;

% traces of the scene parameter estimates accross EM iterations
s_x_trc = zeros(num_sources,1);
s_y_trc = zeros(num_sources,1);
s_b_trc = zeros(num_sources,1);

% log probability for all source position
lnP = log(prob_fn(X(:),Y(:)));
lnP = lnP(:,mode_counts >0);

count = 0;
while ( ~isempty(find(s_x-s_x_trc(:,end), 1)) || ~isempty(find(s_y-s_y_trc(:,end), 1))) && count<=n_em_max  
    
    % increment count
    count = count + 1;
    
    s_x_trc(:,count) = s_x;
    s_y_trc(:,count) = s_y;
    s_b_trc(:,count) = s_b;
    
    % EXPECTATION STEP
    
    % get modal PMFs for each of the source position estimates 
    p_j_est = prob_fn(s_x,s_y); 
      
    % measurement weights
    p_c = s_b .* p_j_est(:,mode_counts>0);
    T = mode_counts(mode_counts>0) .* p_c ./ sum(p_c,1);

    % get Q
    Q = sum(pagemtimes(reshape(T,[size(T,1),1,size(T,2)]), reshape(lnP,[1,size(lnP)])),3);
    
    
    % reshape Q for 2D    
    Q_2D = zeros([size(X),num_sources]);
    for i = 1:num_sources
        Q_2D(:,:,i) = reshape(Q(i,:),size(X));
    end 
    
    %{
    % debugging
    % visualization of objective function for MLE estimators
    peaked_Q = zeros([size(X),num_sources]);
    for i = 1:num_sources
        ki =  find(Q(i,:) == max(Q(i,:))); % max indices
        peaked_Qi = reshape(Q(i,:),size(X));
        peaked_Qi(ki) = peaked_Qi(ki) + 1e4;
        peaked_Q(:,:,i) = peaked_Qi; 
    end 
    %}
    
    % MAXIMIZATION STEP
    %s_b = sum(T,2)/size(T,2); % update source weights
    
    % get max likelihood location indices
    % the following ensures distinct source locations are chosen
    ind_sxy = getMLESourceIndices(Q_2D);
    
    %dx = X(1,2) - X(1,1);
    %dy = Y(2,1) - Y(1,1);
    
    % assign source positions 
    s_x = X(ind_sxy);% + dx*randn(numel(ind_sxy),1); %(plus a bit of noise to avoid degeneracies) 
    s_y = Y(ind_sxy);% + dy*randn(numel(ind_sxy),1); %(plus a bit of noise to avoid degeneracies) 
    
end


end
