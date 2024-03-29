function Visualize_GramSchmidt(nj,mj,X,Y,rl,phi_nm)
    
    nn = unique(nj);
    mm = unique(mj);

    
    figure
    t = tiledlayout(numel(nn),numel(mm));
    %t = tiledlayout(4,1);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    
    
    for j = 1:numel(nj)
        
        nexttile(sub2ind([numel(nn),numel(mm)],nj(j)+1,mj(j)+1))
        %subplot(numel(nn),numel(mm),j)
        
        imagesc([min(X(:)),max(X(:))]/rl,[min(Y(:)),max(Y(:))]/rl,abs(phi_nm(:,:,j)).^2)
        title(['(',num2str(nj(j)),',',num2str(mj(j)),')'])
        axis square

        %xticks([1,ceil(size(phi_nm,2)/2),size(phi_nm,2)]);
        %xticklabels([-1/2,0,1/2])
        xlabel('$x/\sigma$','interpreter','latex')

        %yticks([1,ceil(size(phi_nm,1)/2),size(phi_nm,1)]);
        %yticklabels([-1/2,0,1/2])
        ylabel('$y/\sigma$','interpreter','latex')
    end
    
end
