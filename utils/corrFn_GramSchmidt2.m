function Gamma_nm = corrFn_GramSchmidt2(x,y,Kx,Ky,d2k,GS_basis_mom,A_tot)


   Gamma_nm = 2*pi/sqrt(A_tot) * conj(Basis_GramSchmidt2(x,y,Kx,Ky,d2k,GS_basis_mom));

end
