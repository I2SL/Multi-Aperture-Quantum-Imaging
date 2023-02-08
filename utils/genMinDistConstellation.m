function xy_coords = genMinDistConstellation(n,min_sep)
    % Generates the coordinates of a constellation of n_src point-sources 
    % located in a disk of radius .5 wherein the minimum separation accross
    % all pairs of sources is exactly min_sep. The returned coordinates are 
    % shifted so that the centroid of the constellation (center-of-mass) 
    % lies at the origin (0,0).
    
    assert( 1<n && n <= 20);
    
    D = 1;
    R = D/2;
       
    % packing fractions for circle-within-a-circle up to 20 circles
    % https://en.wikipedia.org/wiki/Circle_packing_in_a_circle
    circ_pack_frac = [1,0.5000,0.6466,0.6864,0.6854,0.6666,0.7777,0.7328,0.6895,0.6878,0.7148,0.7392,0.7245,0.7474,0.7339,0.7512,0.7403,0.7609,0.8034,0.7623];
    
    % no possible samples if the area of the area from the optimal packing fraction
    % is less than the area that the sources may be.
    assert(n*(min_sep/2).^2 < circ_pack_frac(n) *  R^2); 
    
    % generate random sources in the cartesian coordinate system
    xy_coords = rand(n,2) - R;
    
    % rescale so that closest points are min_sep away from each other
    xy_coords = xy_coords * min_sep / min(pdist(xy_coords));
    
    % discard and regenerate points that are outside of the FOV or violated
    % the min_sep distance
    d = pdist(xy_coords);
    p = sum(xy_coords.^2,2);
    dd = d >= min_sep;          % points that pass min sep criterion
    pp = p <= R^2;               % points that pass FOV criterion
    
    % rejection sampling
    while ~all(dd) || ~all(pp)
        
        % resample points that violate the criterion
        xy_coords = rand(n,2) - .5;
        xy_coords = xy_coords * min_sep / min(pdist(xy_coords));
        
        % get new point criterion
        d = pdist(xy_coords);
        p = sum(xy_coords.^2,2);
        dd = d >= min_sep;
        pp = p <= R^2;
    end
    
    % recenter to centroid (note, this may cause some of the points to fall
    % outside of th FOV but it is unlikely),
    %xy_coords = xy_coords - mean(xy_coords,1);
    
end