function p = ModalProb_GramSchmidt_mom(xy_coords,Kx,Ky,d2k,GS_basis_mom,A_tot)
    % Computes the modal probability of detecting a photon in the
    % Gram-Schmidt basis for sources positioned at [x,y].
    
    % correlation function
    correlation_fn = corrFn_GramSchmidt_mom(xy_coords,Kx,Ky,d2k,GS_basis_mom,A_tot); 
    
    %p = abs(correlation_fn).^2 * dx * dy;
    p = abs(correlation_fn).^2;
end