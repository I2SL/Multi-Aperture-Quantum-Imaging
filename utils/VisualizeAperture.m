function R_eff = VisualizeAperture(aperture)

    addpath('utils/circles_v1.1.1')
    
    % multi-aperture parameters
    ap_num = size(aperture,1);
    aper_coords = aperture(:,1:2);
    aper_rads = aperture(:,3); 

    % get the effective aperture diameter
    if ap_num>1
        B = squareform(pdist(aper_coords));             % sub-aperture pairwise centroid distances
        D = aper_rads + aper_rads';                     % sub-aperture pairwise radius sums
        assert(all(triu(B,1) >= triu(D,1),'all'));      % check if any apertures overlap

        % set the effective aperture diameter to the minimum enclosing circle diameter
        cm_coords = aper_coords - mean(aper_coords,1);                          % get centered coordinates (with respect to centroid of centroids -- still need to prove with this works)
        tangent_pts = cm_coords + cm_coords.*aper_rads./vecnorm(cm_coords,2,2); % candidate tangent points where the enclosing circle might touch
        [c_x,c_y,R_eff] = SmallestEnclosingCircle(tangent_pts(:,1)',tangent_pts(:,2)'); % effective aperture radius
        D_eff = 2*R_eff;
    else
        R_eff = aper_rads(1);
        D_eff = 2*R_eff;                         % set the effective aperture diameter to that of the input aperture
    end
    
    
    % plot origin
    scatter(0,0,10,'black','filled')
    % plot the effective aperture
    %circles(c_x,c_y,R_eff,'edgecolor','k','facecolor',[0.0078 0.5765 0.5255],'facealpha',.5)
    circles(c_x,c_y,R_eff,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',.5)
    % plot the sub-aperture
    circles(aperture(:,1),aperture(:,2),aperture(:,3),'facecolor','black')
    
    % figure features
    xlim([-R_eff,R_eff])
    ylim([-R_eff,R_eff])
    axis 'equal'
    %{
    legend({'','Effective Aperture','Compound Aperture'})
    xlabel('K_x')
    ylabel('K_y')

    %}
     
end