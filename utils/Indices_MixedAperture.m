function [nj,mj,vj] = zernikeIndices(n_max,n_apertures)
    % computes the local aperture Zernike mode indices for a system of
    % n apertures.
    
    nj = [];
    mj = [];

    for n = 0:n_max
        for m = -n:2:n
            nj = [nj, n];
            mj = [mj, m];
        end
    end

    n_modes = numel(nj);

    % radial, azimuthal, and aperture index lists
    [vj,jj] = ndgrid(1:n_apertures,1:n_modes);
    nj = nj(jj(:));
    mj = mj(jj(:));
    vj = vj(:)';
    
end
