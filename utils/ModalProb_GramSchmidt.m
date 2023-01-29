function p = ModalProb_GramSchmidt(x,y,X,Y,GS_basis_pos,A_tot)
    % Computes the modal probability of detecting a photon in the
    % Gram-Schmidt basis for sources positioned at [x,y].
    %
    % x,y - source locations
    % X,Y - image plane meshgrids
    % GS_modes  - a matrix stack representing the GS modes over X,Y
    % A_tot - total area of the aperture
    
    % correlation function
    correlation_fn = corrFn_GramSchmidt(x,y,X,Y,GS_basis_pos,A_tot); 
    
    p = abs(correlation_fn).^2;
    %p = [p, max(0,1-sum(p,2))];         % add bucket mode due to truncation
    %p = p ./ sum(p,2);                  % normalize
end
