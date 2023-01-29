function CFI_r_nm = CFI_r_GramSchmidt2(alpha_vec,Kx,Ky,d2k,GS_basis_mom,A_tot,s_b)
    
    % 2-source coordinates
    alpha2_vec = [alpha_vec;-alpha_vec];

    % probability of GS mode nm 
    P_nm = s_b.' * ModalProb_GramSchmidt2(alpha2_vec,Kx,Ky,d2k,GS_basis_mom,A_tot);
    
    % partial derivative of the probability with respect to the
    % half-separation
    dr_P_nm = 2 * s_b.' * real( conj(corrFn_GramSchmidt2(alpha2_vec,Kx,Ky,d2k,GS_basis_mom,A_tot) .* dr_corrFn_GramSchmidt2(alpha2_vec,Kx,Ky,d2k,GS_basis_mom,A_tot)) );
    
    
    % 2-point source separation CFI by mode
    CFI_r_nm = (dr_P_nm).^2 ./ P_nm;
    
end