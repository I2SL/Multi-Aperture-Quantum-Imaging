function dr_Gamma_nm_xy = dr_corrFn_GramSchmidt(xy_coords,Kx,Ky,d2k,GS_basis_mom,A_tot)
    % computes the derivative of the correlation function for the GS basis
    % with respect to the radial coordinate r.
    dr_Gamma_nm_xy = 2*pi/sqrt(A_tot) * conj(dr_GramSchmidt(xy_coords,Kx,Ky,d2k,GS_basis_mom));
end