% Generates a box of locations for possible human motion
% Armspan: https://censusatschool.ca/data-results/2006-07/average-arm-span/

space = rand(10000,3);
space(:,1) = space(:,1)*1.8 - 0.9;
space(:,3) = space(:,3)*1.8 - 0.9;
space(:,2) = space(:,2)*2.5;

idx = (abs(space(:,1)) > 0.2) + (abs(space(:,3)) > 0.2) + (space(:,2) > 1.85);
idx = idx > 0;
gt_data = space(idx,:);

save('body_space.mat', 'gt_data','-v7.3');
%scatter3(space(:,1),space(:,3),space(:,2))