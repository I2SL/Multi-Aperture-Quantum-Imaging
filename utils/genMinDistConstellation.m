function xy = genMinDistConstellation(b, min_sep, align_centroid)
    % Generates the coordinates of a constellation of n_src point-sources 
    % located in a disk of radius .5 wherein the minimum separation accross
    % all pairs of sources is exactly min_sep. 
    % If align_centroid is not 0, the returned coordinates for a  
    % constellation whos centroid (center-of-mass) lies at the origin (0,0).
    % ------------------
    % INPUTS:
    % ------------------
    % b              : nx1 vector containing relative point weights
    % min_sep        : minimum separation distance in constellation
    % align_centroid : require centroid to be at coordinate origin
    % ------------------
    % OUTPUTS
    % ------------------
    % xy             : nx2 matrix containing the xy coordinate pairs of the
    %                  constellation. 
    %                       xy(:,1) is the x coordinate
    %                       xy(:,2) is the y coordinate 
    
    n = numel(b);
    assert((abs(sum(b)- 1) < 1e-16) && all(b>0));
    assert( 1<n && n <= 20);
    
    D = 1;          % Diameter of support disk
    R = D/2;        % Radius of support disk
    
    % packing fractions for circle-within-a-circle up to 20 circles
    % https://en.wikipedia.org/wiki/Circle_packing_in_a_circle
    circ_pack_frac = [1,0.5000,0.6466,0.6864,0.6854,0.6666,0.7777,0.7328,0.6895,0.6878,0.7148,0.7392,0.7245,0.7474,0.7339,0.7512,0.7403,0.7609,0.8034,0.7623];

    % no possible samples if the area of the area from the optimal packing fraction
    % is less than the area that the sources may be.
    assert(n*(min_sep/2).^2 < circ_pack_frac(n) *  R^2); 

    % the point coordinates
    xy = zeros(n,2);
    
    % uniformly sample a point anywhere within the disk  of radius R-min_sep 
    % (ensures that the second point is guaranteed to lie within the disk
    % of radius R)
    p1 = sampleDiskUniform(1,R-min_sep); 
    xy(1,:) = p1;

    % sample a point on the circle of radius min_sep
    [r12_x,r12_y] = pol2cart(rand(1)*2*pi, min_sep); r12 = [r12_x,r12_y];
    
    % place the next point at minimum separation radius from the first
    p2 = p1+r12;
    xy(2,:) = p2;
    
    % if only 2 sources in the constellation then we are done.
    if n == 2
        % recenter the constellation if centroid must be at origin
        if align_centroid
            xy = xy - sum(b.*xy,1);
        end
        return
    end


    % sample the remaining points uniformly from the disk of radius R
    % until reaching a valid constellation (rejection sampling)
    
    isvalid= 0; % indicator for valid constellation
    while ~isvalid  
        % sample the remaining candidate points
        xy(3:n,:) = sampleDiskUniform(n-2,R);

        % recenter the constellation if centroid must be at origin
        if align_centroid
            xy = xy - sum(b.*xy,1);
        end

        % pairwise point distances
        d = pdist(xy);
        
        % a constellation is valid if all points are at least min_sep away
        % from each other and lie within the disk of radius R        
        isvalid = all( (d >= min_sep) & (sum(xy.^2,2) < R^2) );
    end
 
end

function xy = sampleDiskUniform(n,R)
   % generates n samples uniformly over the unit disk
   r = R*sqrt(rand(n,1));
   th = 2*pi*rand(n,1);
   
   [x,y] = pol2cart(th,r);
   xy = [x,y];
end
