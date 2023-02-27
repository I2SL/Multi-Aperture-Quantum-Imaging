function Visualize_MixedAperture(nj,mj,vj,X,Y,rl,U,aperture)
       % A visualization script for seeing the mixed aperture modes.
       % Each aperture index vj is displayed in a new figure page.
 
       num_modes = numel(nj);
       psi_nmv = abs(Basis_MixedAperture([X(:),Y(:)],nj,mj,vj,U,aperture)).^2;
       psi_nmv = reshape(psi_nmv,[size(X),num_modes]);
       
             
       max_n = max(nj);
       sz = [max_n+1,2*max_n+1];
       
       v = 0;
       for j = 1:num_modes
           % trigger a new figure if the aperture index changes
           if vj(j) ~= v
               figure
               v = vj(j);
           end
           
           n = nj(j); m = mj(j);
           k = sub2ind([sz(2),sz(1)],m+max_n+1,n+1);
           subplot(sz(1),sz(2),k)
           imagesc([min(X(:)),max(X(:))]/rl,[min(Y(:)),max(Y(:))]/rl,psi_nmv(:,:,j))
           axis 'square'
           title(['$(',num2str(n),',',num2str(m),')_',num2str(v),'$'],'interpreter','latex')
           xticks([])
           yticks([])
       end
end

%{
function Visualize_MixedAperture(nj,mj,vj,X,Y,rl,U,aper_coords)
       % A visualization script for seeing the mixed aperture modes.
       % Each aperture index vj is displayed in a new figure page.
 
       num_modes = numel(nj);
       psi_nmv = abs(Basis_MixedAperture([X(:),Y(:)],nj,mj,vj,U,aper_coords)).^2;
       psi_nmv = reshape(psi_nmv,[size(X),num_modes]);
       
             
       max_n = max(nj);
       sz = [max_n+1,2*max_n+1];
       
       v = 0;
       for j = 1:num_modes
           % trigger a new figure if the aperture index changes
           if vj(j) ~= v
               figure
               v = vj(j);
           end
           
           n = nj(j); m = mj(j);
           k = sub2ind([sz(2),sz(1)],m+max_n+1,n+1);
           subplot(sz(1),sz(2),k)
           imagesc([min(X(:)),max(X(:))]/rl,[min(Y(:)),max(Y(:))]/rl,psi_nmv(:,:,j))
           axis 'square'
           title(['$(',num2str(n),',',num2str(m),')_',num2str(v),'$'],'interpreter','latex')
           xticks([])
           yticks([])
       end
end
%}