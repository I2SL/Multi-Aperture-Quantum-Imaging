function [nj,mj,vj] = Indices_MixedAperture(n_max,n_apertures)
    % computes the local aperture Zernike mode indices for a system of
    % n apertures.
    
    nj = [];
    mj = [];
    vj = [];
    
    % aperture, radial, and azimuthal index list     
    for v = 1:n_apertures
        for n = 0:n_max
            for m = -n:2:n
                nj = [nj, n];
                mj = [mj, m];
                vj = [vj, v];
            end
        end
    end    
end
