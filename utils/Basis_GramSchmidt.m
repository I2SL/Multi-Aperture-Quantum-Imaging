function psi_nm = Basis_GramSchmidt(x,y,X,Y,GS_basis_pos)
    
    dx = X(1,2) - X(1,1);
    dy = Y(2,1) - Y(1,1);

    n_modes = size(GS_basis_pos,3);
    psi_nm = zeros([numel(x),n_modes]);
    
    for i = 1:n_modes
        GS_i = GS_basis_pos(:,:,i);
        %psi_nm(:,i) = interp2(X,Y,GS_i,x,y) * dx * dy;
        psi_nm(:,i) = interp2(X,Y,GS_i,x,y);
    end
end