function CFI_r_nm_v = CFI_r_MixedAperture(alpha_vec,n,m,v,U,aper_coords,A_tot,b_s)

    
    % Probability of mixed mode (n,m,v)
    P_k_nu = sum(b_s .* abs(corrFn_MixedAperture(alpha_vec,n,m,v,U,aper_coords,A_tot)).^2);
    
    % partial derivative of the probability with respect to the separation
    
    dr_P_k_nu_s = real(   conj(corrFn_MixedAperture(alpha_vec,n,m,v,U,aper_coords,A_tot)) .* dr_corrFn_MixedAperture(alpha_vec,n,m,v,U,aper_coords,A_tot) +...
                        conj(corrFn_MixedAperture(-alpha_vec,n,m,v,U,aper_coords,A_tot)).* dr_corrFn_MixedAperture(-alpha_vec,n,m,v,U,aper_coords,A_tot)...
                     );
            
    % weight by the source brightnesses
    %dr_P_k_nu = sum(b_s .* dr_P_k_nu_s);
    dr_P_k_nu = dr_P_k_nu_s; 
    
    % 2-point source separation CFI
    CFI_r_nm_v = (dr_P_k_nu).^2 ./ P_k_nu;

end