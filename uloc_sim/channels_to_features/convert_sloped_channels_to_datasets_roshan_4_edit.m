% This code creates a dataset for the compensation metric
% We want to create a large dataset and divide it into training and test
% sets.
clear
for index = 1:6
    disp(index)
    load(['/media/extHDD01/roshan_data/channels_pix2pix_disc_',num2str(index),'.mat']) 
    clear channels
    for idx=1:10:31
        disp(idx);
        %% Define Constants

        input_grid_size=0.2;  %Grid size for input Gaussian
        output_grid_size=0.2;   %Grid size for output Gaussian
        output_sigma = 0.25; %Sigma for Gaussian output image
        %ap_orient=[1,1,2,2];% 1 if ap is along x, 2 otherwise
        theta_vals = -pi/2:0.01:pi/2; %Theta values to consider for multipath profiles
        d_vals = -15:0.1:15; %D vals to consider for distance profiles
    %     n_points=length(channels);
        n_lambda=length(cur_model.lambda);
        n_ap=length(ap);
        n_ant=length(ap{1});

        n_features=n_ap; % One distance profile per antenna, one angle profile per ap
        max_x =max(cur_model.walls(:,1));
        max_y =max(cur_model.walls(:,2));
        min_x =min(cur_model.walls(:,1));
        min_y =min(cur_model.walls(:,2));
        d1 = -6.5:input_grid_size:19.5;
        d2 = -2.5:input_grid_size:7.5;
        n_points=100;%*6*6*6*6; % We want equal correct and incorrect points. 
                                                % We multiply by 6, because 4 APs=>6 pairs
                                                % Also, we want to have features in
                                                % both orders.
        choices = nchoosek(1:4,2);
        opt.freq = 3e8./(cur_model.lambda);
        opt.ant_sep = abs(ap{1}(2,1)-ap{1}(1,1));
        features_temp = zeros(n_points,6*6*6*6,n_ap,length(d2),length(d1));
        labels_temp = zeros(n_points,6*6*6*6);
        %% Compensated Examples
        % We want to apply random +-10 m offsets and see if the network can catch
        tau_space = permn([0:5]*cur_model.Ts/5,4);
        tau = linspace(0,cur_model.Ts,6);
        n=0;
        parfor channel_idx=1:100
            [features_temp(channel_idx,:,:,:,:),labels_temp(channel_idx,:)] = compensate_with_given_offsets(channels_offset(100*(idx-1)+channel_idx,:,:,:),offset(100*(idx-1)+channel_idx,:),cur_model,n_ant,opt,d1,d2,n_ap,theta_vals,d_vals,ap);
            if (mod(channel_idx,10)==0)
                disp(channel_idx);
            end
        end
        features = reshape(features_temp,100*6*6*6*6,n_ap,length(d2),length(d1));
        
        
        
%         parfor tau_all = 1:size(tau_space,1)
%             
%             features_temp(channel_idx,:,:,:) = compensate_with_given_offsets(channels_offset(100*idx+channel_idx,:,:,:),tau_space,n_ant,opt);
% %             compensation(:,:,1) = repmat(exp(-1j*opt.freq.'*2*pi/3e8*tau_space(tau_all,1)),1,n_ant);
% %             compensation(:,:,2) = repmat(exp(-1j*opt.freq.'*2*pi/3e8*tau_space(tau_all,1)),1,n_ant);
% %             compensation(:,:,3) = repmat(exp(-1j*opt.freq.'*2*pi/3e8*tau_space(tau_all,1)),1,n_ant);
% %             compensation(:,:,4) = repmat(exp(-1j*opt.freq.'*2*pi/3e8*tau_space(tau_all,1)),1,n_ant);
% %             parfor channel_idx=1:100
% % 
% %                 comp_channels = squeeze(channels_offset(100*idx+channel_idx,:,:,:)).*compensation;
% %                 features_temp(channel_idx,:,:,:) = generate_features_abs(comp_channels,ap,theta_vals,d_vals,d1,d2,opt);
% % 
% %             end
%             n_features_temp = size(features_temp,1);
%             features(n+1:n+n_features_temp,:,:,:) = features_temp;
%             n = n+n_features_temp;
%             disp(n);
%         end
        feature_mean =zeros(1,size(features,2));
        feature_std = zeros(1,size(features,2));
        for k=1:length(feature_mean)
            dset = features(:,k,:,:);
            feature_mean(k) = mean(dset(:));
            feature_std(k) = std(dset(:));
        end
        features = (features-repmat(feature_mean,size(features,1),1,size(features,3),size(features,4)))./repmat(feature_std,size(features,1),1,size(features,3),size(features,4));
        %     Ts_offset = median(diff(tau));
        %     offsets_discrete = Ts_offset.*(round(offsets./Ts_offset));
        labels = reshape(labels_temp,100*6*6*6*6,1);
        save(['/media/extHDD01/roshan_data/dataset_disc_net1_',num2str(index),'_4_',num2str(idx),'.mat'],'features','labels','offset','-v7.3');
    end
end