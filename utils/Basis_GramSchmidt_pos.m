function psi_nm = Basis_GramSchmidt_pos(xy_query,X,Y,GS_basis_pos)
    % computes the GS basis function at the query points (xq,yq) within the
    % discritezation of the image plane defined by X,Y via 2D
    % interpolation.
       
    xq = xy_query(:,1);
    yq = xy_query(:,2);
    
    n_modes = size(GS_basis_pos,3);
    psi_nm = zeros([numel(xq),n_modes]);  
    
    
    for i = 1:n_modes
        GS_i = GS_basis_pos(:,:,i);
        psi_nm(:,i) = interp2(X,Y,GS_i,xq,yq);
    end
    
end