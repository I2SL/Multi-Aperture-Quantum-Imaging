function coords_xy = Polygon(n,rotation,type,param)
    % generates the vertex coordinates of a regular polygon 
    % with n vertices equally spaced around a disk. The coordinates can
    % either lie on a circule of a given radius or be separated 
    % by a given distance,
    %
    % n         : the number of vertices
    % rotation  : the (counter-clockwise) rotation angle of the polygon [radians]
    % type      : ['radius','separation']
    % param     : radius or separation based on type
    
    assert(n>=0, 'Invalid nuber of polygon vertices')
    
    % polygon of order 1 is just a point
    if n == 1
        coords_xy = [0,0];
        return
    end
    
    
    % polygon vertex coordinates around a unit disk of radius 1
    r = ones([n,1]);
    th = 2*pi*(0:n-1)'/n + pi/2 + rotation;
       
    
    switch type
        case 'radius'
            % Radius of vertices 
            R = param;
 
        case 'separation'
            % Radius to make adjacent vertex spacing equal to S 
            S = param;
            del_th = th(2)-th(1);
            
            if n == 2
                R = S/2;
            elseif n > 2
                R = S/2 / sin(del_th/2);
            end
            
    end
    
    r = R*r;
                
    % source positions on a disk of radius R
    [x,y] = pol2cart(th,r); 
    coords_xy = [x,y];
    



end