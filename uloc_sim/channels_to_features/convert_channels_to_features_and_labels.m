clear;

%Load Channel file and initialize variables
load('datasets/channels_1.mat'); 
input_grid_size=1;  %Grid size for input Gaussian
output_grid_size=1;   %Grid size for output Gaussian
output_sigma = 2; %Sigma for Gaussian output image
ap_orient=[1,1,2,2];% 1 if ap is along x, 2 otherwise
theta_vals = [-pi/2:0.01:pi/2]; %Theta values to consider for multipath profiles
d_vals = 0:0.1:25; %D vals to consider for distance profiles
labels_discrete = ceil(labels./output_grid_size);
n_xlabels = length(unique(labels_discrete(:,1)));
n_ylabels =length(unique(labels_discrete(:,2)));
n_points=10000;%size(channels,1);
base_idx=0;
n_lambda=size(channels,2);
n_ap=size(channels,3);
n_ant=size(channels,4);

n_features=n_ap+n_ap*n_ant; % One distance profile per antenna, one angle profile per ap
max_x =max(cur_model.walls(:,1));
max_y =max(cur_model.walls(:,2));
d1 = -max_x:input_grid_size:2*max_x;
d2=-max_y:input_grid_size:2*max_y;
features = zeros(n_points,n_features,length(d2),length(d1));
%labels_gaussian = zeros(size(labels,1),n_xlabels*n_ylabels);
map_X = repmat((1:n_xlabels)'*output_grid_size,1,n_ylabels);
map_Y = repmat((1:n_ylabels)*output_grid_size,n_xlabels,1);

for i=1:n_points
    % Get features
    cur_channels = squeeze(channels(i+base_idx,:,:,:));
    feature_idx=1;
    for j=1:n_ap
        P_angle = compute_multipath_profile(squeeze(cur_channels(1,j,:)),ap{j}(:,ap_orient(j)),cur_model.lambda(1),theta_vals);
        P_out = convert_multipath_to_2d(P_angle,theta_vals,d1,d2,ap{j});
        features(i,feature_idx,:,:)=abs(P_out);
        feature_idx = feature_idx + 1;
        for k=1:n_ant
            P_d = compute_distance_profile(squeeze(cur_channels(:,j,k)),cur_model.lambda,2,d_vals);            
            P_out=convert_distance_to_2d(P_d,d_vals,d1,d2,ap{j}(k,:));
            features(i,feature_idx,:,:)=abs(P_out);
            feature_idx = feature_idx + 1;
        end
    end
    % Get labels
    
%     d = (map_X-labels(i,1)).^2+(map_Y-labels(i,2)).^2;
%     cur_gaussian = exp(-d/output_sigma/output_sigma);%*1/sqrt(2*pi)/output_sigma;
%     cur_gaussian = cur_gaussian.';
%     labels_gaussian(i,:)=cur_gaussian(:);
    if(mod(i,100)==0)
        disp(i);
    end
end
labels_gaussian_2d=[];%get_gaussian_labels(labels(base_idx+1:base_idx+n_points,:),output_grid_size,output_sigma);
    
save('dataset_distance_aoa_large_field_1.mat','features','labels','labels_gaussian_2d','cur_model','ap','-v7.3');


