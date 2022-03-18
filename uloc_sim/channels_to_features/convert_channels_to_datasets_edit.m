clear;
for j=2:7
    disp(j);
    str = ['/media/user1/easystore/datasets/pix2pix/channels/channels_pix2pix_real9x5_',num2str(j),'.mat'];

    %Load Channel file and initialize variables
    load(str); 
%     load('datasets/channels_20_4_rand.mat');
    input_grid_size=0.2;  %Grid size for input Gaussian
    output_grid_size=0.2;   %Grid size for output Gaussian
    output_sigma = 0.5; %Sigma for Gaussian output image
    %ap_orient=[1,1,2,2];% 1 if ap is along x, 2 otherwise
    theta_vals = -pi/2:0.01:pi/2; %Theta values to consider for multipath profiles
    d_vals = -5:0.1:20; %D vals to consider for distance profiles
    n_points=size(channels_w_offset,1);
    n_lambda=length(cur_model.lambda);
    n_ap=length(ap);
    n_ant=length(ap{1});

    n_features=n_ap; % One distance profile per antenna, one angle profile per ap
    max_x =max(cur_model.walls(:,1));
    max_y =max(cur_model.walls(:,2));
    min_x =min(cur_model.walls(:,1));
    min_y =min(cur_model.walls(:,2));
%     d1 = min_x:input_grid_size:max_x;
%     d2 = min_y:input_grid_size:max_y;
    d1 = -4.5:input_grid_size:13.5;
    d2 = -2.5:input_grid_size:7.5;
%     labels_gaussian = zeros(size(labels,1),n_xlabels*n_ylabels);
%     map_X = repmat((1:n_xlabels)'*output_grid_size,1,n_ylabels);
%     map_Y = repmat((1:n_ylabels)*output_grid_size,n_xlabels,1);

%     channels_rel= zeros(length(cur_model.lambda),n_ap,n_ant,n_ap-1);
    opt.freq = 3e8./(cur_model.lambda);
%     opt2.freq = 3e8./(cur_model2.lambda);
    opt.ant_sep = 0.0259;
    opt.threshold = 0.1;
%     opt2.ant_sep = abs(ap{1}(1,2)-ap{1}(2,2));
    
    features_wo_offset = zeros(n_points,n_ap,length(d2),length(d1));
    features_w_offset = zeros(n_points,n_ap,length(d2),length(d1));
    parfor i=1:n_points
        % Get features
        features_wo_offset(i,:,:,:) =generate_features_abs(squeeze(channels_wo_offset(i,:,:,:)),ap,theta_vals,d_vals,d1,d2,opt);
        features_w_offset(i,:,:,:) = generate_features_abs(squeeze(channels_w_offset(i,:,:,:)),ap,theta_vals,d_vals,d1,d2,opt);
        
        features_wo_offset_filt(i,:,:,:) =generate_features_abs_filt(squeeze(channels_wo_offset(i,:,:,:)),ap,theta_vals,d_vals,d1,d2,opt);
        features_w_offset_filt(i,:,:,:) = generate_features_abs_filt(squeeze(channels_w_offset(i,:,:,:)),ap,theta_vals,d_vals,d1,d2,opt);
        % Get labels
    %         d = (map_X-labels(i,1)).^2+(map_Y-labels(i,2)).^2;
    %         cur_gaussian = exp(-d/output_sigma/output_sigma)*1/sqrt(2*pi)/output_sigma;
    %         labels_gaussian(i,:)=cur_gaussian(:);
        if(mod(i,1000)==0)
            disp(i);
        end    
        
%         figure(1)
%         for k=1:4
%             subplot(4,2,k), imagesc(d1,d2,squeeze(features_wo_offset(i,k,:,:))), axis xy, hold on...
%                 ,scatter(labels(i,1),labels(i,2),'r','filled'), hold off;
%             subplot(4,2,4+k), imagesc(d1,d2,squeeze(features_wo_offset_music(i,k,:,:))), axis xy, hold on...
%                 ,scatter(labels(i,1),labels(i,2),'r','filled'), hold off;
%         end
%         waitforbuttonpress;
    end
    %%
    label_std = 1/sqrt(2*pi);%std(labels_gaussian_all(:));
    feature_mean =zeros(1,size(features_wo_offset,2));
    feature_std = zeros(1,size(features_wo_offset,2));
    for k=1:length(feature_mean)
        dset = features_wo_offset(:,k,:,:);
        feature_mean(k) = mean(dset(:));
        feature_std(k) = std(dset(:));
    end
    features_wo_offset = (features_wo_offset-repmat(feature_mean,size(features_wo_offset,1),1,size(features_wo_offset,3),size(features_wo_offset,4)))./repmat(feature_std,size(features_wo_offset,1),1,size(features_wo_offset,3),size(features_wo_offset,4));
%     labels_gaussian_2d=get_gaussian_labels1(labels_discrete,output_grid_size,2);
    
    feature7_mean =zeros(1,size(features_w_offset,2));
    feature7_std = zeros(1,size(features_w_offset,2));
    for k=1:length(feature7_mean)
        dset1 = features_w_offset(:,k,:,:);
        feature7_mean(k) = mean(dset1(:));
        feature7_std(k) = std(dset1(:));
    end
    features_w_offset = (features_w_offset-repmat(feature7_mean,size(features_w_offset,1),1,size(features_w_offset,3),size(features_w_offset,4)))./repmat(feature7_std,size(features_w_offset,1),1,size(features_w_offset,3),size(features_w_offset,4));
    %%
    feature_mean =zeros(1,size(features_wo_offset_filt,2));
    feature_std = zeros(1,size(features_wo_offset_filt,2));
    for k=1:length(feature_mean)
        dset = features_wo_offset_filt(:,k,:,:);
        feature_mean(k) = mean(dset(:));
        feature_std(k) = std(dset(:));
    end
    features_wo_offset_filt = (features_wo_offset_filt-repmat(feature_mean,size(features_wo_offset_filt,1),1,size(features_wo_offset_filt,3),size(features_wo_offset_filt,4)))./repmat(feature_std,size(features_wo_offset_filt,1),1,size(features_wo_offset_filt,3),size(features_wo_offset_filt,4));
%     labels_gaussian_2d=get_gaussian_labels1(labels_discrete,output_grid_size,2);
    
    feature7_mean =zeros(1,size(features_w_offset_filt,2));
    feature7_std = zeros(1,size(features_w_offset_filt,2));
    for k=1:length(feature7_mean)
        dset1 = features_w_offset_filt(:,k,:,:);
        feature7_mean(k) = mean(dset1(:));
        feature7_std(k) = std(dset1(:));
    end
    features_w_offset_filt = (features_w_offset_filt-repmat(feature7_mean,size(features_w_offset_filt,1),1,size(features_w_offset_filt,3),size(features_w_offset_filt,4)))./repmat(feature7_std,size(features_w_offset_filt,1),1,size(features_w_offset_filt,3),size(features_w_offset_filt,4));
    
    %%
    
    %labels_gaussian_2d=get_gaussian_labels2(labels_discrete,output_grid_size,1,length(d1),length(d2));
    labels_gaussian_2d=get_gaussian_labels_negative(labels,output_grid_size,output_sigma,d1,d2);
    stri = ['/media/user1/easystore/datasets/pix2pix/features/datasets_pix2pix_real9x5_2x_',num2str(j),'.mat'];

    save(stri,'features_w_offset','features_wo_offset','labels','labels_gaussian_2d','cur_model','ap','-v7.3');
    
    stri = ['/media/user1/easystore/datasets/pix2pix/features/datasets_pix2pix_real9x5_2xfilt_',num2str(j),'.mat'];

    save(stri,'features_w_offset_filt','features_wo_offset_filt','labels','labels_gaussian_2d','cur_model','ap','-v7.3');
end