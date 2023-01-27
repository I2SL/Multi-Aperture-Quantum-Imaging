function Gamma_nm_v = corrFn_MixedAperture(xy_coords,n,m,v,U,aper_coords,A_tot)
       % correlation function for the Mixed Aperture Basis
       Gamma_nm_v =  2*pi/sqrt(A_tot) * Basis_MixedAperture(xy_coords,n,m,v,U,aper_coords);
end