clear;

for j=1:6
    str = ['datasets/channel20_10_',num2str(j),'.mat'];
    %Load Channel file and initialize variables
    load(str); 
    
    input_grid_size=0.2;  %Grid size for input Gaussian
    output_grid_size=0.2;   %Grid size for output Gaussian
    output_sigma = 1; %Sigma for Gaussian output image
    ap_orient=[1,1,2,2];% 1 if ap is along x, 2 otherwise
    theta_vals = -pi/2:0.01:pi/2; %Theta values to consider for multipath profiles
    d_vals = -15:0.1:15; %D vals to consider for distance profiles
    n_points=10000;
    n_lambda=length(cur_model.lambda);
    n_ap=length(ap);
    n_ant=length(ap{1});

    n_features=n_ap*n_ap; % One distance profile per antenna, one angle profile per ap
    max_x =max(cur_model.walls(:,1));
    max_y =max(cur_model.walls(:,2));
    d1 = 0:input_grid_size:max_x;
    d2 = 0:input_grid_size:max_y;
    model7 = load_model_edit(7);
%     labels_gaussian = zeros(size(labels,1),n_xlabels*n_ylabels);
%     map_X = repmat((1:n_xlabels)'*output_grid_size,1,n_ylabels);
%     map_Y = repmat((1:n_ylabels)*output_grid_size,n_xlabels,1);

    channels_rel= zeros(length(cur_model.lambda),n_ap,n_ant,n_ap-1);
    opt.freq = 3e8./(cur_model.lambda);
    opt.ant_sep = abs(ap{1}(2,1)-ap{1}(1,1));
    
    features = zeros(n_points,n_ap*n_ap,length(d2),length(d1));
    features7 = zeros(n_points,n_ap*n_ap,length(d2),length(d1));
    i=0;
    parfor i=1:n_points
        % Get features
        features(i,:,:,:) = generate_features_new(squeeze(channels(i,:,:,:)),ap,theta_vals,d_vals,d1,d2,opt,cur_model);
        features7(i,:,:,:) = generate_features_new(squeeze(channels(i,:,:,:)),ap,theta_vals,d_vals,d1,d2,opt,model7);
        % Get labels
    %         d = (map_X-labels(i,1)).^2+(map_Y-labels(i,2)).^2;
    %         cur_gaussian = exp(-d/output_sigma/output_sigma)*1/sqrt(2*pi)/output_sigma;
    %         labels_gaussian(i,:)=cur_gaussian(:);
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
%     labels_gaussian_2d=get_gaussian_labels1(labels_discrete,output_grid_size,2);
    
    feature7_mean =zeros(1,size(features7,2));
    feature7_std = zeros(1,size(features7,2));
    for k=1:length(feature7_mean)
        dset1 = features7(:,k,:,:);
        feature7_mean(k) = mean(dset1(:));
        feature7_std(k) = std(dset1(:));
    end
    features7 = (features7-repmat(feature7_mean,size(features7,1),1,size(features7,3),size(features7,4)))./repmat(feature7_std,size(features7,1),1,size(features7,3),size(features7,4));
    labels_gaussian_2d=get_gaussian_labels1(labels_discrete,output_grid_size,1);
    
    stri = ['datasets/dataset_GAN20_',num2str(j),'.mat'];
    save(stri,'features','features7','labels_discrete','labels_gaussian_2d','cur_model','ap','-v7.3');
    clear
end