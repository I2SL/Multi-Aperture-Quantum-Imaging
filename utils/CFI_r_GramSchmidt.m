function CFI_r_nm = CFI_r_GramSchmidt(xy_coords,X,Y,GS_basis_pos,A_tot,s_b)
    x=xy_coords(:,1); y = xy_coords(:,2);
    
    % GS correlation function
    p_nm = sum(s_b .* ModalProb_GramSchmidt(x,y,X,Y,GS_basis_pos,A_tot),1);
    
    % Take the numerical derivative of the probability wrt to the radial coordinate
    [theta,~] = cart2pol(x,y);
    dr = sqrt(2)*abs(X(1,2)-X(1,1));
    dx = dr*cos(theta);
    dy = dr*sin(theta);
    dr_p_nm_s = (ModalProb_GramSchmidt(x+dx,y+dy,X,Y,GS_basis_pos,A_tot) - p_nm) / dr;
    dr_p_nm = sum(s_b.*dr_p_nm_s,1);
    
    % the CFI
    CFI_r_nm = (dr_p_nm).^2 ./ p_nm;
    
end
