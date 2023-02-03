function [s_x,s_y] = MLESourceCoords(X,Y,Q)
    % Returns the MLE source coordinates for the supplied Q function.
    %
    % X,Y - Coordinate matrices of dimension NxN
    % Q - matrix stack of dimensions NxNxsrc_num

    src_num = size(Q,3);    % number of sources
    
    H = cell(1,src_num); % a cell array of candidate coordinates for each source
    I = cell(1,src_num);
    
    for i = 1:src_num
        s_ix = X(Q(:,:,i) == max(Q(:,:,i),[],'all'));
        s_iy = Y(Q(:,:,i) == max(Q(:,:,i),[],'all'));
        H{i} = [s_ix,s_iy];   % candidate positions for source i
        I{i} = 1:numel(s_ix); % number of candidates positions for source i
    end
    
    % get all combinations of source coordinate indices
    D = I;
    [D{:}] = ndgrid(I{:});
    Z = cell2mat(cellfun(@(m)m(:),D,'uni',0)); % a matrix where each row is a unique index set into the coordinates in H
    
    
    % make all constellations
    spread = zeros(size(Z,1),1);
    constellations = zeros(src_num,2,numel(spread));
    for j = 1:numel(spread)
        for k = 1:src_num
            hk = H{k};
            constellations(k,:,j) = hk(Z(j,k),:);
        end
        spread(j) = sum(1./pdist(constellations(:,:,j)));  % spread of the constellation (metric modelled as potential energy of a system of point charges)
    end
    
    % select constellations witht the largest spreads
    % (the potential energy of the charge configuration is a minimum).
    c = (spread == min(spread));    % candidate constellation indices
    s_x = constellations(:,1,c);
    s_y = constellations(:,2,c);

end