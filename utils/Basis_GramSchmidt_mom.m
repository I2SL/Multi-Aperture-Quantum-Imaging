function GS_basis_pos_xy = Basis_GramSchmidt_mom(xy_coords,Kx,Ky,d2k,GS_basis_mom)
    % Computes the Gram-Schmidt basis in position space at the sample locations x,y
    % using an explicit continuous inverse Fourier Transform.
    x = xy_coords(:,1);
    y = xy_coords(:,2);
    
    GS_basis_pos_xy = ctsIFT_2D(x,y,Kx,Ky,d2k,GS_basis_mom);

end
