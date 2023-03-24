function [s_x,s_y] = MLESourceCoords(X,Y,Q,src_coords,method)
% Returns the MLE source coordinates for the supplied Q function.
    %
    % X,Y - Coordinate matrices of dimension NxN
    % Q - matrix stack of dimensions NxNxsrc_num
    % src_coords - ground truth constellation coordinates
    % method - ['Heuristic','MinError']

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
    

    % make all candidate MLE constellations
    nc = size(Z,1); % number of candidate constellations
    constellations = zeros(src_num,2,nc);
    for j = 1:nc
        for k = 1:src_num
            hk = H{k};
            constellations(k,:,j) = hk(Z(j,k),:);
        end
    end
    
    % choose among the candidate MLE constellations using the provided method
    switch method
        case 'Heuristic'
            % uses a heuristic to choose among degenerate solutions
            % (this is effectively folding in some prior information about the scene generation mechanism)
            
            % model the sources as electrostatic charges and compute the potential energy of the configuration
            charge_PE = zeros(size(Z,1),1); 
            for j = 1:nc
                charge_PE(j) = sum(1./pdist(constellations(:,:,j)));  % potential energy of a system of point charges
            end
            
            % constellations with the greates "spread" are the minimum electrostatic energy configurations
            spread = 1./charge_PE;
            c = (spread == max(spread));    % candidate constellation indices
            sel_i = c(randi(1:numel(c)));   % randomly select one of the candidates that maximize the heuristic
            s_x = constellations(:,1,sel_i);
            s_y = constellations(:,2,sel_i);
            
        case 'MinError'
            % select based on minimum error solution (this is cheating a bit)
            loc_err = zeros(1,nc);
            for i = 1:nc   
                loc_err(i) = LocalizationError(src_coords,constellations(:,:,i));
            end
            [~, sel_i] = min(loc_err);  
            s_x = constellations(:,1,sel_i);
            s_y = constellations(:,2,sel_i);
            
            %{
            % visualize
            figure
            scatter(src_coords(:,1),src_coords(:,2),'filled','black')
            hold on
            for i = 1:nc 
                scatter(constellations(:,1,i),constellations(:,2,i),'filled');
            end
            %}
    end
    
    
    %{
    num_constellations = size(constellations,3);
    f =  factor(num_constellations);
    half_point = floor(numel(f)/2);
    figure
    tiledlayout(prod(f(1:half_point)),prod(f(half_point+1:end)))
    for c = 1:size(constellations,3)
        nexttile
        scatter(constellations(:,1,c),constellations(:,2,c),'filled')
        axis equal
        title('Energy: ',num2str(spread(c)));
    end
    %}

end