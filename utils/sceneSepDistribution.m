num_src = 5;
num_trials = 100000;
D = [];
b = ones(num_src,1)/num_src;
min_sep = 1/5;
centroid_aligned = 1;

for t = 1:num_trials
    scene = genMinDistConstellation(b, min_sep, centroid_aligned);
    d = pdist(scene);
    D = [D, d];    
    
    %{
    figure
    scatter(scene(:,1),scene(:,2),'filled','black')
    xlim([-.5,.5])
    ylim([-.5,.5])
    xlabel('x [rl]')
    ylabel('y [rl]')
    title('1/8 rl Min Separation')
    %}
    
end

histogram(D)
title({'Random Root Chain-Link Scene Generation','5 Sources','100K Scenes'})
ylabel('Counts')
xlabel('Source Separation [rl]')