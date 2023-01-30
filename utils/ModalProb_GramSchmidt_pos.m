function p = ModalProb_GramSchmidt_pos(xy_coords,X,Y,GS_basis_pos,A_tot)
    % Computes the modal probability of detecting a photon in the
    % Gram-Schmidt basis for sources positioned at [x,y].
    %
    % x,y - source locations
    % X,Y - image plane meshgrids
    % GS_basis_pos  - a matrix stack representing the GS modes over X,Y
    % A_tot - total area of the aperture
    
    % correlation function
    correlation_fn = corrFn_GramSchmidt_pos(xy_coords,X,Y,GS_basis_pos,A_tot); 
    
    p = abs(correlation_fn).^2;
end
