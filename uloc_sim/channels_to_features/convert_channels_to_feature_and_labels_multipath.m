clear;

%Load Channel file and initialize variables
load('datasets/channels_1.mat'); 
input_grid_size=0.5;  %Grid size for input Gaussian
output_grid_size=2;   %Grid size for output Gaussian
output_sigma = 1; %Sigma for Gaussian output image
ap_orient=[1,1,2,2];% 1 if ap is along x, 2 otherwise
theta_vals = [-pi/2:0.05:pi/2]; %Theta values to consider for multipath profiles
d_vals = -25:1:25; %D vals to consider for distance profiles
labels_discrete = ceil(labels./output_grid_size);
n_xlabels = length(unique(labels_discrete(:,1)));
n_ylabels =length(unique(labels_discrete(:,2)));
n_points=1;%size(channels,1);
n_lambda=size(channels,2);
n_ap=size(channels,3);
n_ant=size(channels,4);

n_features=n_ap*n_ap; % One distance profile per antenna, one angle profile per ap
max_x =max(cur_model.walls(:,1));
max_y =max(cur_model.walls(:,2));
d1 = 0:input_grid_size:max_x;
d2=0:input_grid_size:max_y;
features = zeros(n_points,n_features,length(d2),length(d1));
labels_gaussian = zeros(size(labels,1),n_xlabels*n_ylabels);
map_X = repmat((1:n_xlabels)'*output_grid_size,1,n_ylabels);
map_Y = repmat((1:n_ylabels)*output_grid_size,n_xlabels,1);

channels_rel= zeros(length(cur_model.lambda),n_ap,n_ant,n_ap-1);
opt.freq = 3e8./(cur_model.lambda);
opt.ant_sep = abs(ap{1}(2,1)-ap{1}(1,1));

for i=1:n_points
    % Get features
    cur_channels = squeeze(channels(i,:,:,:));
    feature_idx=1;
    for j=1:n_ap
        for k=1:n_ant
            ind = 1:1:n_ap;
            ind(j) = [];
            for l = 1:n_ap-1
                channels_rel(:,j,k,l) = cur_channels(:,j,k).*conj(cur_channels(:,ind(l),1));
                
            end
        end
        
        for l = 1:n_ap-1
            
            Prel =  compute_multipath_profile2d_fast(squeeze(channels_rel(:,j,:,l)),theta_vals,d_vals,opt);
            Prel_out = convert_relative_spotfi_to_2d_edit(Prel,ap{j},ap{ind(l)},theta_vals,d_vals,d1,d2);
            features(i,feature_idx,:,:) = db(abs(Prel_out));
            feature_idx = feature_idx+1;
            
        end
        P = compute_multipath_profile2d_fast(squeeze(cur_channels(:,j,:)),theta_vals,d_vals,opt);
        P_out = convert_spotfi_to_2d(P,theta_vals,d_vals,d1,d2,ap{j});
        
        features(i,feature_idx,:,:) = db(abs(P_out));
        feature_idx = feature_idx+1;
        
    end
    % Get labels
    
    d = (map_X-labels(i,1)).^2+(map_Y-labels(i,2)).^2;
    cur_gaussian = exp(-d/output_sigma/output_sigma)*1/sqrt(2*pi)/output_sigma;
    labels_gaussian(i,:)=cur_gaussian(:);
    if(mod(i,100)==0)
        disp(i);
    end
end
    
save('dataset_d1_d2_large.mat','features','labels','labels_gaussian','cur_model','ap','-v7.3');


