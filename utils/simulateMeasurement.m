function [mode_count,varargout] = simulateMeasurement(mu_pho,p,varargin)
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
    
    if isPoiss
        mode_count = poissrnd(mu_pho*p);
    else
        mode_count = round(mu_pho*p);
    end
    
        
    % add dark-current photons
    dc = poissrnd(dark_lambda, [1,numel(p)]);
    mode_count = mode_count + dc;
    
    % expand mode counts into measurement
    if nargout == 2
        varargout{1} = repelem(1:numel(p),mode_count);
    end
end

