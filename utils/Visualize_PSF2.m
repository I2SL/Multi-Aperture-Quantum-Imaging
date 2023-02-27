function Visualize_PSF2(X,Y,rl,aperture)

    % absolute square of the PSF
    PSF2 = reshape(abs(MultiAperturePSF([X(:),Y(:)],aperture)).^2,size(X));
        
    % Visualize the psf
    figure

    imagesc([min(X(:)),max(X(:))]/rl,[min(Y(:)),max(Y(:))]/rl,PSF2);
    title('Multi-Aperture |PSF|^2')
    xlabel('x [rl]')
    ylabel('y [rl]')
    axis 'square'
end