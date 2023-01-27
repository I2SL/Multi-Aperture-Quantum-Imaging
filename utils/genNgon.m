function coords_xy = genNgon(n,rand_rotation)
    % generates the vertex coordinates of polygon with n vertices equally spaced
    % around a disk. The vertex separation is defined to be 1.
    if n == 1
        coords_xy = [0,0];
        return
    end
    
    
    r = ones([n,1])/2;
    th = 2*pi*(0:n-1)'/n + pi/2;
    del_th = th(2)-th(1);
    
    if rand_rotation
        th = th + 2*pi*rand(1);
    end
    
    % source positions on a unit disk (radius = 1)
    [x,y] = pol2cart(th,r); 
    coords_xy = [x,y];
    
    % make vertex separation 1 by rescaling the radius    
    if n > 2
        coords_xy = coords_xy / sin(del_th/2);
    end
end