clear
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration used: 3 APs, with switched 8 antenna arrays each
% Wi-Fi channel used: 149, freq: 5745 MHz, Bandwidth: 20 MHz
% Number of subcarriers available: 52, total subcarriers: 64
% Subcarrier order: [-21, -7, 7, 21, -26:-22, -20:-8,-6:-1, 1:6, 8:20, 22:26]
% AP Locations are in ap_pos_20April.mat
% The cell array ap_pos has 3 cells, each containing the position of the 8 antennas in x,y format
% Client locations are in final_locs_20_April.mat in the variable final_locs in x,y format
% If both x and y location is zero, ignore this measurement
% Channels are in channels_20April.mat
% The variable h has all the channels. h is a N_measurement times 3 cell array
% In each cell, there is 8 by 52 array of channels, the 8 rows correspond to 8 antennas and 52 columns correspond to 52 subcarriers in the order specified above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data
% load('channels_20April.mat')
% load('final_locs_20_April.mat')
% load('ap_pos_20April.mat')
[subbands,subband_order] = sort([-21, -7, 7, 21, -26:-22, -20:-8,-6:-1, 1:6, 8:20, 22:26]); % channels to considered
opt.freq = (5745 + (subbands).*(20/64))*10^6; % all the subcarrier frequencies
% % freq2 = (2060 + (-64:63).*(50/128))*10^6; % all the subcarrier frequencies
opt.lambda = 3e8./opt.freq;
% % lambda2 = 3e8./freq2;
% n_h = length(h);
% channel = zeros(440,52,8,3);
% count=0;
% indexes = [];
% temp_pos=final_locs;
% for i=1:n_h
%     if(~(final_locs(i,1)==0 && final_locs(i,2)==0))
%         count = count+1;
%         indexes = [indexes,i];
%         parfor k=1:3
%             channel(count,:,:,k) = h{i,k}';
%         end
%         cli_pos(count,:) = final_locs(i,:);
%     end
% end
% channel = channel(:,subband_order,:,:);
% 
% ap = ap_pos
load('real_data_3ant_20Apr.mat');
%% Feature Calculation
opt.lambda = 3e8./opt.freq;
n_lambda=length(opt.lambda);
n_points = size(channel,1);
n_ap=length(ap_pos);
n_ant=length(ap_pos{1});
theta_vals = -pi/2:0.01:pi/2; %Theta values to consider for multipath profiles
% d_vals = -8:0.1:8;
% d_vals1 = 0:0.1:5;
opt.ap_orient = [1,2,2];
% feature_idx=1;
% opt.freq = [freq1(channels),freq2(channels)];
opt.ant_sep = abs(ap_pos{1}(2,1)-ap_pos{1}(1,1));
d1=-2:0.1:2;
d2=-2:0.1:3;
% % features = zeros(n_points,n_ap*(n_ap-1),length(d2),length(d1));
features = zeros(n_points,n_ap,length(d2),length(d1));
% i=0;
parfor i=1:n_points
%     % Get features
    features(i,:,:,:) = generate_features_abs_aoa(squeeze(channel(i,:,:,:)),ap_pos,theta_vals,d1,d2,opt);
    if(mod(i,1000)==0)
        disp(i);
    end    
end
label_std = 1/sqrt(2*pi);%std(labels_gaussian_all(:));
feature_mean =zeros(1,size(features,2));
feature_std = zeros(1,size(features,2));
for k=1:length(feature_mean)
    dset = features(:,k,:,:);
    feature_mean(k) = mean(dset(:));
    feature_std(k) = std(dset(:));
end
features = (features-repmat(feature_mean,size(features,1),1,size(features,3),size(features,4)))./repmat(feature_std,size(features,1),1,size(features,3),size(features,4));

labels_gaussian_2d=get_gaussian_labels_negative(cli_pos,0.1,0.25,d1,d2);
labels_discrete = round(cli_pos*10)/10;

% save('datasets/dataset2_real27Feb_ext.mat','features','labels_gaussian_2d','labels_discrete','-v7.3');