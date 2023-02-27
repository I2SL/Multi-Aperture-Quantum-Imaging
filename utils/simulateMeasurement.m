function [measurement,mode_count] = simulateMeasurement(mu_pho,p,isPoiss)
    % mu_pho - mean photon number 
    % p - modal probability distribution (PMF)
    % isPoiss - trigger for sampling the the number of collected photons
    %           from a poisson distribution

    % set the number of photons collected in the measurement
    if isPoiss
        N = poissrnd(mu_pho); 
    else
        N = mu_pho;
    end
    
    % randomly assign modal bin to each photon according to PMF
    modes = 1:numel(p);
    measurement = datasample(modes,N,'Weights',p);
    
    % count photons in modal bins
    [gc,gs] = groupcounts(measurement');
    mode_count = zeros(1,numel(p));
    mode_count(gs) = gc;
end

