function [measurement,mode_count] = simulateMeasurement(mu_pho,p,varargin)
    % mu_pho - mean photon number 
    % p - modal probability distribution (PMF)
    % 
    % ----------------
    % Optional Inputs:
    % ----------------
    % isPoiss - trigger for sampling the the number of collected photons
    %           from a poisson distribution
    % dark_lambda - poisson rate of dark current at each photon counter
    
    
    %%%%%%%%%%%%%%%%% Parser %%%%%%%%%%%%%%%%%%%%%%%%%%
    default_isPoiss = 1;
    default_dark_lambda = 0;
    
    P = inputParser;
    
    addRequired(P,'mu_pho');
    addRequired(P,'p');
    
    addOptional(P,'isPoiss',default_isPoiss);
    addOptional(P,'dark_lambda',default_dark_lambda);
    
    parse(P, mu_pho, p, varargin{:});
    
    isPoiss = P.Results.isPoiss;
    dark_lambda = P.Results.dark_lambda;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % set the number of photons collected in the measurement
    if isPoiss
        N = poissrnd(mu_pho); 
    else
        N = mu_pho;
    end
    
    % randomly assign modal bin to each photon according to PMF
    modes = 1:numel(p);
    measurement = datasample(modes,N,'Weights',p);
    
    % add dark-current photons
    dc = poissrnd(dark_lambda, [1,numel(p)]);
    measurement = [measurement,repelem(1:numel(p),dc)];
    
    % count photons in modal bins
    [gc,gs] = groupcounts(measurement');
    mode_count = zeros(1,numel(p));
    mode_count(gs) = gc;
end

