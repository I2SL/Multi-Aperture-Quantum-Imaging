function CFI_r_nm = CFI_r_GramSchmidt(alpha_vec,Kx,Ky,d2k,GS_basis_mom,A_tot,s_b)
    
    % 2-source coordinates
    alpha2_vec = [alpha_vec;-alpha_vec];
    
    
    % correlation function
    corr_fn = corrFn_GramSchmidt_mom(alpha2_vec,Kx,Ky,d2k,GS_basis_mom,A_tot);
    
    % radial derivative of correlation function
    dr_corr_fn = dr_corrFn_GramSchmidt(alpha2_vec,Kx,Ky,d2k,GS_basis_mom,A_tot);

    % probability of GS mode nm 
    P_nm = s_b.' * abs(corr_fn).^2;
    
    % radial derivative of the probability 
    dr_P_nm = 2 * s_b.' * real( conj(corr_fn) .* dr_corr_fn );
    
    % 2-point source separation CFI by mode
    CFI_r_nm = (dr_P_nm).^2 ./ P_nm;
    
end