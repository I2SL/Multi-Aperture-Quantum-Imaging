function Visualize_PSF2(X,Y,rl,aper_coords)
    % Visualize the psf    
    figure
    PSF2 = reshape(abs(MultiAperturePSF([X(:),Y(:)],aper_coords)).^2,size(X));
    imagesc([min(X(:)),max(X(:))]/rl,[min(Y(:)),max(Y(:))]/rl,PSF2);
    title('Multi-Aperture |PSF|^2')
    xlabel('x [rl]')
    ylabel('y [rl]')
    axis 'square'
end