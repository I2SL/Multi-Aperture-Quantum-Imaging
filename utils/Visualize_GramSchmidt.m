function VisualizeModes_GramSchmidt(nj,mj,phi_nm)
    
    nn = unique(nj);
    mm = unique(mj);

    
    figure
    for j = 1:numel(nj)
        
        subplot(numel(nn),numel(mm),j)
        
        imagesc(abs(phi_nm(:,:,j)).^2)
        title(['(',num2str(nj(j)),',',num2str(mj(j)),')'])
        axis square

        xticks([1,ceil(size(phi_nm,2)/2),size(phi_nm,2)]);
        xticklabels([-1/2,0,1/2])
        xlabel('rl')

        yticks([1,ceil(size(phi_nm,1)/2),size(phi_nm,1)]);
        yticklabels([-1/2,0,1/2])
        ylabel('rl')
    end
    
end
