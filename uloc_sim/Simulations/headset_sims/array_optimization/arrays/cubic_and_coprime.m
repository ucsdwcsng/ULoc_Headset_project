%% Generate and Evaluate Cubic and Coprime Arrays
addpath('../../')
clearvars
c = 3e8;
opt.freq = 4.4928e9; opt.freq = 6489.6e6;
opt.lambda = c./opt.freq;

% Generate Cubic Array
ant_sep = opt.lambda *3/8;
dimensions = [0.165,0.090,0.070]; % width height x depth in meters of Oculus Quest 2
n = floor(dimensions./ant_sep);
n = [11, 7, 7];
x = ant_sep*[0:n(1)-1];
y = ant_sep*[0:n(2)-1];
z = ant_sep*[0:n(3)-1];
[X,Y,Z] = meshgrid(x,y,z);
x1 = [X(:),Y(:),Z(:)];
x1 = x1 - mean(x1);
figure(1)
scatter3(x1(:,1),x1(:,3),x1(:,2))
title('Cubic Array'); view(-85,30);
axis equal; xlabel('Width'); ylabel('Depth'); zlabel('Height');
S1 = gen_auto_corr_steering(deg2rad(-180:5:180), deg2rad(-90:5:90), opt, x1.');

% Generate Coprime Array
% Closest Coprime array is slightly larger than the above, needing antenna
% size [10, 6, 6]
n = [10, 6, 6];
M = [2,5]; N = [2,3]; J = [2,3];
Mv{1} = (0:M(1))*ant_sep*M(2); Mv{2}=(0:M(2))*ant_sep*M(1);
Nv{1} = (0:N(1))*ant_sep*N(2); Nv{2}=(0:N(2))*ant_sep*N(1);
Jv{1} = (0:J(1))*ant_sep*J(2); Jv{2}=(0:J(2))*ant_sep*J(1);
[X,Y,Z] = meshgrid(Mv{1},Nv{1},Jv{1});
x2 = [X(:),Y(:),Z(:)];
[X,Y,Z] = meshgrid(Mv{2},Nv{2},Jv{2});
x2 = [x2;X(:),Y(:),Z(:)];
x2 = unique(x2,'rows');
x2 = x2-mean(x2);
figure(2)
subplot(2,2,1)
scatter3(x2(:,1),x2(:,3),x2(:,2))
title('Coprime Array Top'); view(0,90);
axis equal; xlabel('Width'); ylabel('Depth'); zlabel('Height');
subplot(2,2,2)
scatter3(x2(:,1),x2(:,3),x2(:,2))
title('Coprime Array Front'); view(0,0);
axis equal; xlabel('Width'); ylabel('Depth'); zlabel('Height');
subplot(2,2,3)
scatter3(x2(:,1),x2(:,3),x2(:,2))
title('Coprime Array Iso'); view(45,45);
axis equal; xlabel('Width'); ylabel('Depth'); zlabel('Height');
subplot(2,2,4)
scatter3(x2(:,1),x2(:,3),x2(:,2))
title('Coprime Array Side'); view(90,0);
axis equal; xlabel('Width'); ylabel('Depth'); zlabel('Height');
S2 = gen_auto_corr_steering(deg2rad(-180:5:180), deg2rad(-90:5:90), opt, x2.');

% Autocorrelations of each
figure(3)
subplot(2,2,1)
imagesc(abs(S1))
title('Cubic Autocorrelation: 539 Antennas')
axis tight equal

subplot(2,2,2)
imagesc(abs(S2))
title('Coprime Autocorrelation: 115 Antennas')
axis tight equal
subplot(2,1,2)
imagesc(S2-S1)
title('Difference Heatmap: Coprime-Correlation')
axis tight equal
colorbar


