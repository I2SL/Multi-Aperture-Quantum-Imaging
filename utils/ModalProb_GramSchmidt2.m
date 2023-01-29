function p = ModalProb_GramSchmidt2(xy_coords,Kx,Ky,d2k,GS_basis_mom,A_tot)
    % Computes the modal probability of detecting a photon in the
    % Gram-Schmidt basis for sources positioned at [x,y].
    
    % correlation function
    correlation_fn = corrFn_GramSchmidt2(xy_coords,Kx,Ky,d2k,GS_basis_mom,A_tot); 
    
    %p = abs(correlation_fn).^2 * dx * dy;
    p = abs(correlation_fn).^2;
    %p = [p, max(0,1-sum(p,2))];         % add bucket mode due to truncation
    %p = p ./ sum(p,2);                  % normalize
end