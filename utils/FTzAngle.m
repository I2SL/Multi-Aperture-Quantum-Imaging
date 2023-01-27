function v = FTzAngle(theta,m)
    % Computes the angular function of the Fourier Transformed Zernikes
    
    v = zeros(numel(theta),numel(m));
    
    % angular function
    c = cos(abs(m).*theta);
    s = sin(abs(m).*theta);
    
    v(:,m>0) = sqrt(2) * c(:,m>0);
    v(:,m==0)= 1;
    v(:,m<0) = sqrt(2) * s(:,m<0);
end

