function GS_modes_pos_xy = GramSchmidtBasis2(x,y,j_nm,Kx,Ky,aperture,GS_modes_mom)
    % Manually evaluates the GS basis in position space at the sample location x,y
    % to avoid the implicit 

    d2k = (Kx(1,2)-Kx(1,1))*(Ky(2,1)-Ky(1,1));
 

    GS_modes_pos_xy = 1/(2*pi) * sum( GS_modes_mom(:,:,j_nm) .* aperture .* ...
                                      exp(1i* ( Kx.*reshape(x,[1,1,numel(x)]) + ...
                                                Ky.*reshape(y,[1,1,numel(y)]) ) ) ,[1,2])...
                                            .* d2k;

end
