function dr_Gamma_nm_v = dr_corrFn_MixedAperture(xy_coords,n,m,v,U,aper_coords,A_tot)
% derivative of Gamma wrt to source half-separation
[theta,r] = cart2pol(xy_coords(:,1),xy_coords(:,2));
dr_Gamma_nm_v =  2*pi/sqrt(A_tot) * conj(dr_FTZernike(r,theta,n,m) .* phaseFn(xy_coords,v,U,aper_coords)...
                                     + FTZernike(r,theta,n,m) .* dr_phaseFn(xy_coords,v,U,aper_coords));    
end
