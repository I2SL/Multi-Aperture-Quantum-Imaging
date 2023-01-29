function CFI_r_nmv = CFI_r_MixedAperture(alpha_vec,n,m,v,U,aper_coords,A_tot,s_b)
    % 2-source coordinates
    alpha2_vec = [alpha_vec;-alpha_vec];
    
    % Probability of mixed mode (n,m,v)
    P_nmv = s_b.' * ModalProb_MixedAperture(alpha2_vec,n,m,v,U,aper_coords,A_tot);
    
    % partial derivative of the probability with respect to the 
    % half-separation                
    dr_P_nmv = 2 * s_b.' * real( conj(corrFn_MixedAperture(alpha2_vec,n,m,v,U,aper_coords,A_tot)) .* dr_corrFn_MixedAperture(alpha2_vec,n,m,v,U,aper_coords,A_tot) );

    % 2-point source separation CFI by mode
    CFI_r_nmv = (dr_P_nmv).^2 ./ P_nmv;

end