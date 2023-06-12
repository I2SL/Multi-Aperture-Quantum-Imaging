D = 30;         % multi-aperture effective diameter [length]
R = D/2;    % multi-aperture effective radius   [length]
d = 3;              % sub-aperture diameter             [length]
r = d/2;            % sub-apeture radius                [length]

mono = [0,0,R];
plus9 = [PlusAperture(9,R-r),r*ones(9,1)];
ring9 = [Polygon(9,0,'radius',R-r),r*ones(9,1)];
golay9 = [Golay9(R-r),r*ones(9,1)];    

apertures = {mono,plus9,ring9,golay9};

fig = figure;
fig.RendererMode = 'manual';
fig.Renderer = 'painters';
t = tiledlayout(1,numel(apertures));
%t.TileSpacing = 'compact';
%t.Padding = 'compact';

for a = 1:numel(apertures)
    nexttile
    VisualizeAperture(apertures{a})
    axis off
end

