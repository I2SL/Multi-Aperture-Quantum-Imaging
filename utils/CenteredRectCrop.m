function I_out =  CenteredRectCrop(I,rect_r,rect_c)
    % crops the matrix I to a rectangle of dimensions [rect_r, rect_c]
    % that is centered on the matrix.
    
    row_max = size(I,1);
    col_max = size(I,2); 
    
    center = [ceil(row_max/2),ceil(col_max/2)];
    delta = [floor(rect_r/2),floor(rect_c/2)];
    
    I_out = I( (center(1)-delta(1)) : (center(1)+delta(1)) ,...
               (center(2)-delta(2)) : (center(2)+delta(2)));
end