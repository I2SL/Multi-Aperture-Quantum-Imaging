function [measurement,mode_count] = simulateMeasurement(n_pho,p)
    N = poissrnd(n_pho); % number of photons collected
    
    % randomly assign modal bin to each photon according to PMF
    modes = 1:numel(p);
    measurement = datasample(modes,N,'Weights',p);
    
    % count photons in modal bins
    [gc,gs] = groupcounts(measurement');
    mode_count = zeros(1,numel(p));
    mode_count(gs) = gc;
end

