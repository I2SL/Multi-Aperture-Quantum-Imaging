function ShowScene(scene,est_scene)
% Plots the constellation and the estimate (if available).
% scene         - nx3 matrix characterizing the scene. The columns are organized as [x [rl], y [rl], brightness]
% est_scene     - nx3 matrix characterizing the estimate. The columns are organized as [x [rl], y [rl], brightness]

% number of sources
num_src = size(scene,1);  

% plot
figure
hold on
scatter(scene(:,1),scene(:,2),36*num_src*scene(:,3),'filled','black')
if ~isempty(est_scene)
    scatter(est_scene(:,1),est_scene(:,2),36*num_src*est_scene(:,3),'+','filled','red')
end
hold off
title([num2str(num_src),'-Source Configuration'])
xlabel('x [rl]')
ylabel('y [rl]')
xlim([-.5,.5])
ylim([-.5,.5])
xticks(linspace(-.5,5,5))
yticks(linspace(-.5,5,5))

xticklabels({'-1/2','-1/4','0','1/4','1/2'})
yticklabels({'-1/2','-1/4','0','1/4','1/2'})
if ~isempty(est_scene)
    legend({'scene','estimate'})
else
    legend({'scene'})
end

end