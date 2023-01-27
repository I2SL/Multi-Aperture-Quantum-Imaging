function VisualizeModes_MixedAperture(nj,mj,vj,X,Y,U,aper_coords)
    
 
       num_modes = numel(nj);
       psi_nmv = abs(Basis_MixedAperture([5*X(:),5*Y(:)],nj,mj,vj,U,aper_coords)).^2;
       psi_nmv = reshape(psi_nmv,[size(X),num_modes]);
       
       pages = unique(vj);
       modes = max(nj)+1;
       [n,inj,in] = unique(nj);
       [m,imj,im] = unique(mj);
       
       
       
       figs = cell(numel(pages),1);
       for p = pages
           figs{p} = figure(p);           
       end
       
       
       sz = [max(in),max(im)];
       for j = 1:num_modes
           k = sub2ind(sz,in(j),im(j));
           %set(0,'CurrentFigure',figs{p}) 
           subplot(sz(1),sz(2),k,'Parent', figs{p})
           imagesc(psi_nmv(:,:,j))
           axis 'square'
           title(['(',num2str(n(in(j))),',',num2str(m(im(j))),')'])
           
       end

end