function Gamma_nm = corrFn_GramSchmidt2(xy_coords,Kx,Ky,d2k,GS_basis_mom,A_tot)
    % correlation function for the Gram-Schmidt basis
    Gamma_nm = 2*pi/sqrt(A_tot) * conj(Basis_GramSchmidt2(xy_coords,Kx,Ky,d2k,GS_basis_mom));
end
