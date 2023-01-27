function dr_Gamma_nm = dr_corrFn_GramSchmidt(x,y,theta,GS_basis_mom,Kx,Ky,A_tot)
    % Derivative of GS correlation function with respect to radial coordinate
    d2k = (Kx(1,2)-Kx(1,1))*(Ky(2,1)-Ky(1,1));
    dr_Gamma_nm = 1/A_tot * sum( GS_basis_mom .* (1i*cos(theta)*Kx + 1i*sin(theta)*Ky) .* ...
                                exp(1i * ( Kx .* reshape(x,[1,1,numel(x)]) +...
                                           Ky .* reshape(y,[1,1,numel(y)]) ) ), [1,2] ) * d2k;

end