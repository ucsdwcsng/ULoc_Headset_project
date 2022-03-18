%% Evaluate metrics from headset_sim
% Given saved mat files for "metrics" and "array" from headset_sim.m:
clearvars
dir = ['metrics/'];

m{1} = load([dir,'axial_3-17.mat']);
m{2} = load([dir,'zz_3-17.mat']);
m{3} = load([dir,'opt_3-17.mat']);

labels{1} = 'Axial';
labels{2} = 'ZigZag';
labels{3} = 'Optimized';

%% AoA Errors
for ii=1:length(m)
    figure(1)
    subplot(2,2,1)
    hold on
    ecdf(abs(m{ii}.metrics.AoA_err(:,1)))
    labels_prct1{ii} = [labels{ii},', 50: ',num2str(prctile(abs(m{ii}.metrics.AoA_err(:,1)),50)), ', 90: ',num2str(prctile(abs(m{ii}.metrics.AoA_err(:,1)),90))];
    subplot(2,2,2)
    hold on
    ecdf(abs(m{ii}.metrics.AoA_err(:,2)))
    labels_prct2{ii} = [labels{ii},', 50: ',num2str(prctile(abs(m{ii}.metrics.AoA_err(:,2)),50)), ', 90: ',num2str(prctile(abs(m{ii}.metrics.AoA_err(:,2)),90))];
    subplot(2,2,3)
    hold on
    ecdf(abs(m{ii}.metrics.close_peak_err(:,1)))
    labels_prct3{ii} = [labels{ii},', 50: ',num2str(prctile(abs(m{ii}.metrics.close_peak_err(:,1)),50)), ', 90: ',num2str(prctile(abs(m{ii}.metrics.close_peak_err(:,1)),90))];
    subplot(2,2,4)
    hold on
    ecdf(abs(m{ii}.metrics.close_peak_err(:,2)))
    labels_prct4{ii} = [labels{ii},', 50: ',num2str(prctile(abs(m{ii}.metrics.close_peak_err(:,2)),50)), ', 90: ',num2str(prctile(abs(m{ii}.metrics.close_peak_err(:,2)),90))];
end
subplot(2,2,1)
grid on; grid minor; title('AoA(Theta) Error'); xlabel('Theta Err (\circ)');ylim([0,1]);xlim([0,60]);
legend(labels_prct1)
subplot(2,2,2)
grid on; grid minor; title('AoA(Phi) Error'); xlabel('Phi Err (\circ)');ylim([0,1]);xlim([0,60]);
legend(labels_prct2)
subplot(2,2,3)
grid on; grid minor; title('Close Peak AoA(Theta) Error'); xlabel('Theta Err (\circ)');ylim([0,1]);xlim([0,60]);
legend(labels_prct3)
subplot(2,2,4)
grid on; grid minor; title('Close Peak AoA(Phi) Error'); xlabel('Phi normed Err (\circ)');ylim([0,1]);xlim([0,60]);
legend(labels_prct4)

%% Num Peaks
figure(2)
h = zeros(length(m),8);
for ii=1:length(m)
    subplot(2,2,1)
    tmp=histogram(m{ii}.metrics.num_peaks);
    pk_idx = unique(m{ii}.metrics.num_peaks);
    pk_idx = min(pk_idx):max(pk_idx);
    try 
        h(ii,pk_idx) = tmp.Values;
        failed = 0;
    catch
        failed = 1;
    end
end
clf
subplot(2,2,1)
hold on
if failed == 1
    for ii=1:length(m)
        histogram(m{ii}.metrics.num_peaks);
    end
else
    bar(h.')
end
grid on; grid minor; grid on; grid minor; title('Number of Peaks per Array'); xlabel('# of peaks'); ylabel('Number Packets');
legend(labels)

%% Eccentricity
figure(2)
for ii=1:length(m)
    subplot(2,2,2)
    hold on
    histogram(m{ii}.metrics.eccentricity);
    %ecc(ii,36-length(tmp.Values):35) = tmp.Values;
end
subplot(2,2,2)
%bar(-0.02:0.03:1,ecc.')
grid on; grid minor; grid on; grid minor; title('Eccentricity per Array'); xlabel('Eccentricity'); ylabel('Number Packets');
legend(labels)

%% Peak Size
figure(2)
subplot(2,2,3)
for ii=1:length(m)
    hold on
    ecdf(m{ii}.metrics.peak_size)
end
grid on; grid minor; grid on; grid minor; title('CDF Peak Size per Array'); xlabel('Peak Size (Pixels)');
legend(labels)

figure(2)
subplot(2,2,4)
for ii=1:length(m)
    hold on
    histogram(m{ii}.metrics.peak_size)
end
grid on; grid minor; grid on; grid minor; title('Hist Peak Size per Array'); xlabel('Peak Size (Pixels)');ylabel('Number Packets');
legend(labels)

%% Surface Plot
theta_step = deg2rad(1);%0.01;
t_vals = -deg2rad(180):theta_step:deg2rad(180); % aoa values for FFT
p_vals = -deg2rad(90):theta_step:deg2rad(90);
[X,Y] = meshgrid(p_vals,t_vals);
ii = 2;
p = zeros(361,181);
pct = p;
for jj = 1:length(m{ii}.metrics.peak_size)
    p(round(m{ii}.metrics.close_peak(jj,1)+181), round(m{ii}.metrics.close_peak(jj,2)+91)) = p(round(m{ii}.metrics.close_peak(jj,1)+181), round(m{ii}.metrics.close_peak(jj,2)+91)) + m{ii}.metrics.peak_size(jj);
    pct(round(m{ii}.metrics.close_peak(jj,1)+181), round(m{ii}.metrics.close_peak(jj,2)+91)) = pct(round(m{ii}.metrics.close_peak(jj,1)+181), round(m{ii}.metrics.close_peak(jj,2)+91)) + 1;
end
pct=pct+(pct==0);
p=p./pct;
figure
s=surf(rad2deg(X),rad2deg(Y),p);title('Peak Size at Theta,Phi');xlabel('Phi Vals');ylabel('Theta Vals');zlabel('Peak Size (pixels)');
s.EdgeColor = 'none';
s=imagesc(rad2deg(t_vals),rad2deg(p_vals),p);title('Peak Size at Theta,Phi');xlabel('Phi Vals');ylabel('Theta Vals');zlabel('Peak Size (pixels)');


%% Arrays
if 1 % Show Arrays
    figure
    X = -0.08*[0 0 0 0 0 1; 1 0 1 1 1 1; 1 0 1 1 1 1; 0 0 0 0 0 1];
    Y = -0.08*[0 0 0 0 1 0; 0 1 0 0 1 1; 0 1 1 1 1 1; 0 0 1 1 1 0];
    Z = -0.08*[0 0 1 0 0 0; 0 0 1 0 0 0; 1 1 1 0 1 1; 1 1 1 0 1 1];
    for ii=1:length(m)
        subplot(1,length(m),ii)
        hold on
        scatter3(m{ii}.array.local_ants{1}(:,1),m{ii}.array.local_ants{1}(:,3),m{ii}.array.local_ants{1}(:,2),...
            30,'filled','r')
        fill3(X,Y,Z,'blue', 'FaceAlpha',0.2)
        view(145,35)
        %view(90,0)
        grid on
        xlabel('X')
        ylabel('Z')
        zlabel('Y')
        title(labels{ii})
        axis equal
    end
end
