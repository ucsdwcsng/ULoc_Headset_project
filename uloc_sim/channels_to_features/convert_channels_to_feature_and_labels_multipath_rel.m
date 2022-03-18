close all;
clear
n = 0;
% j =4;
for j=1:6
    disp(j);
    str = ['datasets/channels_',num2str(j),'.mat'];
    %Load Channel file and initialize variables
    load(str); 
    model = load_model(j);
    input_grid_size=1;  %Grid size for input Gaussian
    output_grid_size=1;   %Grid size for output Gaussian
    output_sigma = 1; %Sigma for Gaussian output image
    ap_orient=[1,1,2,2];% 1 if ap is along x, 2 otherwise
    theta_vals = [-pi/2:0.01:pi/2]; %Theta values to consider for multipath profiles
    d_vals = -15:0.1:15; %D vals to consider for distance profiles
    labels_discrete(n+1:n+20000,:) = ceil(labels(1:20000,:)./output_grid_size);
    n_xlabels = length(unique(labels_discrete(:,1)));
    n_ylabels =length(unique(labels_discrete(:,2)));
    n_points=20000;
    n_lambda=size(channels,2);
    n_ap=size(channels,3);
    n_ant=size(channels,4);

    n_features=n_ap*n_ap; % One distance profile per antenna, one angle profile per ap
    max_x =max(cur_model.walls(:,1));
    max_y =max(cur_model.walls(:,2));
    d1 = 0:input_grid_size:max_x;
    d2=0:input_grid_size:max_y;
    
%     labels_gaussian = zeros(size(labels,1),n_xlabels*n_ylabels);
    map_X = repmat((1:n_xlabels)'*output_grid_size,1,n_ylabels);
    map_Y = repmat((1:n_ylabels)*output_grid_size,n_xlabels,1);

    channels_rel= zeros(length(cur_model.lambda),n_ap,n_ant,n_ap-1);
    opt.freq = 3e8./(cur_model.lambda);
    opt.ant_sep = abs(ap{1}(2,1)-ap{1}(1,1));

    parfor i=1:n_points
        % Get features
        features(n+i,:,:,:) = generate_features_rel_only(pos,ap,theta_vals,d_vals,d1,d2,opt,model);
        % Get labels

%         d = (map_X-labels(i,1)).^2+(map_Y-labels(i,2)).^2;
%         cur_gaussian = exp(-d/output_sigma/output_sigma)*1/sqrt(2*pi)/output_sigma;
%         labels_gaussian(i,:)=cur_gaussian(:);
        if(mod(i,1000)==0)
            disp(i);
        end
    end
    n = n+20000;
end
save('datasets/dataset_d1_d2.mat','features','labels_discrete','cur_model','ap','-v7.3');