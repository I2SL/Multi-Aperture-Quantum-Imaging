function coords_xy = genRandSource(rl,n_sources)
    % generates a random source constellation
    x = rl*(rand(n_sources,1)-0.5);
    y = rl*(rand(n_sources,1)-0.5);
   
    coords_xy = [x,y];
end
