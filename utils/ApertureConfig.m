function [A,Kx,Ky] = ApertureConfig(a_kx,a_ky,ap_dim)
    % returns the normalized multi-aperture function
    % a_kx and a_ky are the coordinates of the sub-apertures in fractional
    % units of the aperture radius (dimensionless).
    
    kx_delta = max(abs(a_kx))+1;
    ky_delta = max(abs(a_ky))+1;
    delta = max(kx_delta,ky_delta);
    
    % aperture plane coordinate grid
    [Kx,Ky] = meshgrid(linspace(-delta,delta,ap_dim));
    dkx = (Kx(1,2)-Kx(1,1)); dky =(Ky(2,1)-Ky(1,1));
    
    circ = @(u,v) (u.^2 + v.^2 < 1);
    
    % construct aperture
    A = zeros(size(Kx));
    for k = 1:numel(a_kx)
        A = A + circ(Kx-a_kx(k),Ky - a_ky(k));
    end
    
    % normalize the aperture
    A = A / sqrt(sum(dkx*dky * abs(A).^2, 'all'));   
end