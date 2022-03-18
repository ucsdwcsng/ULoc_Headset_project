% Inputs: trainpath, testpath, exp_name, 
%         location: Jacobs or Atkinson
% outputs: CDFs
% Plausible path_ids, give them in the same sequence as passed to the
% DeepLoc network:
% train_jacobs_Aug16_1, test_jacobs_Aug16_1
% train_jacobs_Aug16_2, test_jacobs_Aug16_2
% train_jacobs_Aug16_3, test_jacobs_Aug16_3
% train_jacobs_Aug16_4_ref, test_jacobs_Aug16_4_ref
% fov_train_jacobs_July28, fov_test_jacobs_July28,
% non_fov_train_jacobs_July28, non_fov_test_jacobs_July28
% fov_train_jacobs_July28_2, fov_test_jacobs_July28_2,
% non_fov_train_jacobs_July28_2, non_fov_test_jacobs_July28_2
% edit_jacobs_July28, edit_jacobs_July28_2

%The following Inputs are just an example replace them with your Inputs
%trainpath = {'edit_jacobs_July28','fov_train_jacobs_July28_2','non_fov_train_jacobs_July28_2'};
%testpath = {'fov_test_jacobs_July28_2','non_fov_test_jacobs_July28_2'};
%exp_name = 'e0_enc_2dec_train';
clearvars
close all
% trainpath = {'train_atk_vr'};%{, 'train_atk_cart', 'train_jacobs_bot_w_vr', 'train_jacobs_vr'}
% trainpath = {'train_atk_cart'};
trainpath = {'train_jacobs_bot_w_vr'};
% trainpath = {'train_jacobs_vr'};

% testpath = {'test_atk_vr'};%{, 'test_atk_cart', 'test_jacobs_bot_w_vr', 'test_jacobs_vr'}
% testpath = {'test_atk_cart'};
testpath = {'test_jacobs_bot_w_vr'};
% testpath = {'test_jacobs_vr'};

exp_name = 'e304_real2dec_train'; % 'e301_real2dec_train', 'e303_real2dec_train', 'e304_real2dec_train'
env = 'jacobs'; % or 'atk'

%% Define Jacobs hall parameters
ant_sep = 0.0259; % sanatan; ant_sep was not defined for the following lines of code; moved it up here from line 51.
if strcmp(env, 'jacobs')
    ap{1} = [0,2+1.5*ant_sep;...
        0,2+0.5*ant_sep;...
        0,2-0.5*ant_sep;...
        0,2-1.5*ant_sep];
    ap{2} = [16.8,7.6-1.5*ant_sep;...
        16.8,7.6-0.5*ant_sep;...
        16.8,7.6+0.5*ant_sep;...
        16.8,7.6+1.5*ant_sep];
    ap{3} = [6.4+1.5*ant_sep,7.6;...
        6.4+0.5*ant_sep,7.6;...
        6.4-0.5*ant_sep,7.6;...
        6.4-1.5*ant_sep,7.6];
    ap{4} = [12.4-1.5*ant_sep,0;...
        12.4-0.5*ant_sep,0;...
        12.4+0.5*ant_sep,0;...
        12.4+1.5*ant_sep,0];
    max_x =27;
    max_y =12;
    min_x =-9;
    min_y =-4;
elseif strcmp(env, 'atk')
    ap{1} = [0,2.5+1.5*ant_sep;...
        0,2.5+0.5*ant_sep;...
        0,2.5-0.5*ant_sep;...
        0,2.5-1.5*ant_sep];
    ap{2} = [4-1.5*ant_sep,0;...
        4-0.5*ant_sep,0;...
        4+0.5*ant_sep,0;...
        4+1.5*ant_sep,0];
    ap{3} = [8,2.5-1.5*ant_sep;...
        8,2.5-0.5*ant_sep;...
        8,2.5+0.5*ant_sep;...
        8,2.5+1.5*ant_sep];
    ap{4} = [4+1.5*ant_sep,4.5;...
        4+0.5*ant_sep,4.5;...
        4-0.5*ant_sep,4.5;...
        4-1.5*ant_sep,4.5];
    max_x = 12;
    max_y = 7.5;
    min_x = -4;
    min_y = -2.5;
else
    disp('error')
end
scale=0.1;

d1 = min_x:scale:max_x;
d2 = min_y:scale:max_y;
ant_sep = 0.0259;
theta_vals = -pi/2:0.01:pi/2;
d_vals = 0:0.1:30;
subcarrier_indices = [-122:-104,-102:-76,-74:-40,-38:-12,-10:-2,2:10,12:38,40:74,76:102,104:122];% 80MHz
opt.freq = double(5e9 + 5*155*1e6) + subcarrier_indices.*80e6./256;
opt.lambda = 3e8./opt.freq;
opt.ant_sep = ant_sep;
zeroFreq = mean(opt.freq);
lambda_center = 3e8./zeroFreq;
Ts = 1/(80e6);
n_ap = 4;
%% Load appropriate Data
% load Training data
labels_train = [];
for i=1:length(trainpath)
    pathid = trainpath{i};
    switch pathid
        case 'train_atk_vr'
            load('/media/user1/easystore/datasets/aditya_sim/vr_cart_channels/channels_atk_vr.mat','labels');
            n_points_train = round(size(labels,1)*0.70);
            labels_train = [labels_train;labels(n_points_train+1:end,:)];
        case 'train_atk_cart'
            load('/media/user1/easystore/datasets/aditya_sim/vr_cart_channels/channels_atk_cart.mat','labels');
            n_points_train = round(size(labels,1)*0.70);
            labels_train = [labels_train;labels(n_points_train+1:end,:)];
        case 'train_jacobs_bot_w_vr'
            load('/media/user1/easystore/datasets/aditya_sim/vr_cart_channels/channels_jacobs_bot_w_vr.mat','labels');
            n_points_train = round(size(labels,1)*0.85);
            labels_train = [labels_train;labels(n_points_train+1:end,:)];
        case 'train_jacobs_vr'
            load('/media/user1/easystore/datasets/aditya_sim/vr_cart_channels/channels_jacobs_vr.mat','labels');
            n_points_train = round(size(labels,1)*0.85);
            labels_train = [labels_train;labels(n_points_train+1:end,:)];
        otherwise
            error(['The names are non existent. '...
                'Should be on of the following pathnames: '...
                'train_jacobs_Aug16_1 test_jacobs_Aug16_1'...
                'train_jacobs_Aug16_2 test_jacobs_Aug16_2 '...
                'train_jacobs_Aug16_3 test_jacobs_Aug16_3'...
                'train_jacobs_Aug16_4_ref test_jacobs_Aug16_4_ref'...
                'fov_train_jacobs_July28 fov_test_jacobs_July28'...
                'non_fov_train_jacobs_July28 non_fov_test_jacobs_July28'...
                'fov_train_jacobs_July28_2 fov_test_jacobs_July28_2'...
                'non_fov_train_jacobs_July28_2 non_fov_test_jacobs_July28_2'...
                'edit_jacobs_July28 edit_jacobs_July28']);
    end
end
% load Testing data
labels_test = [];
channels_test = [];
features_w_test = [];
features_wo_test = [];
labels_2d_test = [];
for i=1:length(testpath)
    pathid = testpath{i};
    load(['/media/user1/easystore/datasets/aditya_sim/vr_cart_features/dataset_',pathid,'.mat']);
    features_w_test = [features_w_test;features_w_offset];
    features_wo_test = [features_wo_test;features_wo_offset];
    labels_2d_test = [labels_2d_test;labels_gaussian_2d];
    switch pathid
        case 'test_atk_vr'
            load('/media/user1/easystore/datasets/aditya_sim/vr_cart_channels/channels_atk_vr.mat','labels', 'channels_w_offset');
            n_points_train = round(size(labels,1)*0.70);
            labels_test = [labels_test;labels(n_points_train+1:end,:)];
            channels_test = [channels_test;channels_w_offset(n_points_train+1:end,:,:,:)];
        case 'test_atk_cart'
            load('/media/user1/easystore/datasets/aditya_sim/vr_cart_channels/channels_atk_cart.mat','channels_w_offset', 'labels');
            load('/media/user1/easystore/datasets/aditya_sim/vr_cart_channels/channels_atk_vr.mat','labels');
            n_points_train = round(size(labels,1)*0.70);
            labels_test = [labels_test;labels(n_points_train+1:end,:)];
            channels_test = [channels_test;channels_w_offset(n_points_train+1:end,:,:,:)];
        case 'test_jacobs_bot_w_vr'
            load('/media/user1/easystore/datasets/aditya_sim/vr_cart_channels/channels_jacobs_bot_w_vr.mat', 'channels_w_offset', 'labels');
%             load('/media/user1/easystore/datasets/aditya_sim/vr_cart_channels/channels_jacobs_vr.mat','labels');
            n_points_train = round(size(labels,1)*0.85);
            labels_test = [labels_test;labels(n_points_train+1:end,:)];
            channels_test = [channels_test;channels_w_offset(n_points_train+1:end,:,:,:)];
        case 'test_jacobs_vr'
            load('/media/user1/easystore/datasets/aditya_sim/vr_cart_channels/channels_jacobs_vr.mat','labels', 'channels_w_offset');
            n_points_train = round(size(labels,1)*0.85);
            labels_test = [labels_test;labels(n_points_train+1:end,:)];
            channels_test = [channels_test;channels_w_offset(n_points_train+1:end,:,:,:)];
        otherwise
            error(['The names are non existent. '...
                'Should be on of the following pathnames: '...
                'train_jacobs_Aug16_1 test_jacobs_Aug16_1'...
                'train_jacobs_Aug16_2 test_jacobs_Aug16_2 '...
                'train_jacobs_Aug16_3 test_jacobs_Aug16_3'...
                'train_jacobs_Aug16_4_ref test_jacobs_Aug16_4_ref'...
                'fov_train_jacobs_July28 fov_test_jacobs_July28'...
                'non_fov_train_jacobs_July28 non_fov_test_jacobs_July28'...
                'fov_train_jacobs_July28_2 fov_test_jacobs_July28_2'...
                'non_fov_train_jacobs_July28_2 non_fov_test_jacobs_July28_2'...
                'edit_jacobs_July28 edit_jacobs_July28']);
    end
end
clear channels_w_offset features_w_offset features_wo_offset gaussian_labels_2d labels
%% Run SpotFi code
for k=1:size(labels_test,1)
    parfor i=1:n_ap
        P_AoA = abs(compute_multipath_profile2d_fast_edit(squeeze(channels_test(k,:,:,i)),theta_vals,d_vals,opt));
        [~,idx] = max(P_AoA(:)); [aoa_result,d_aoa_result] = ind2sub(size(P_AoA),idx);
        d_aoa_result = d_vals(d_aoa_result); aoa_result=rad2deg(theta_vals(aoa_result));

        aoa_pred_temp(i) = aoa_result;
    end
    aoa_pred(k,:) = aoa_pred_temp;
    if(mod(k,1000)==0)
        disp(k)
    end
end
if strcmp(env, 'jacobs')
    ap_aoa = [0,180,-90,90];
elseif strcmp(env, 'atk')
    ap_aoa = [0,90,180,-90];
end
for j=1:n_ap
    ap_center(j,:) = mean(ap{j},1);
end

[act_error,act_fin_error,act_loc_estimate] = run_Loc_allPoints_edited(ap_center,aoa_pred+ap_aoa,labels_test);
figure(1),cdfplot(act_fin_error), hold on;xlabel('Localization Error(m)'), ylabel('CDF'), hold on;
disp(median(act_fin_error))
%% Run Analyss on the DeepLoc results
decoder = load(['/media/user1/easystore/DeepLoc/results/',exp_name,'/rw_test/decoder.mat'],'outputs');
outputs = double(cell2mat(decoder.outputs.'));
[n_pts,m,n] = size(outputs);

% n_pts = length(decoder.wo_outputs.');
% [n_chan,m,n] = size(decoder.wo_outputs{1,1});
% features_wo_outputs = zeros(n_pts,n_chan,m,n);
% for i=1:n_pts
%     features_wo_outputs(i,:,:,:) = decoder.wo_outputs{1,i};
%     decoder.wo_outputs(1,i) = {zeros(1,1)};
%     if(mod(i,1000)==0)
%         disp(i);
%     end
% end
clear decoder;
figure(1)
%%
parfor i=1:n_pts 
%     temp = squeeze(labels_2d_test(i,:,:)); % labels
%     [~,l_max]=max(temp(:));
    
    map_sel = squeeze(outputs(i,:,:));
    map_xel = map_sel;
    
    temp = map_xel; %predcitions
    [~,p_max]=max(temp(:));
%     [li(i),lj(i)]=ind2sub([m,n],l_max);
    [pi(i),pj(i)]=ind2sub([m,n],p_max);
    
    if (pi(i) < 11)
        pi(i) = 11;
    end
    if (pi(i) > m-11)
        pi(i) = m-11;
    end
    
    if (pj(i) < 11)
        pj(i) = 11;
    end
    if (pj(i) > n-11)
        pj(i) = n-11;
    end
    
    x = d1(pj(i)-10:pj(i)+10);
    y = d2(pi(i)-10:pi(i)+10);
    
    image = squeeze(map_xel(pi(i)-10:pi(i)+10,pj(i)-10:pj(i)+10));
    
    xq = min(x):scale/10:max(x);
    yq = min(y):scale/10:max(y);
    [X,Y] = meshgrid(x,y);
    [Xq,Yq] = meshgrid(xq,yq);
    imageq = interp2(X,Y,image,Xq,Yq,'cubic');
    
    [~,image_max]=max(imageq(:));
    [image_i(i),image_j(i)]=ind2sub(size(imageq),image_max);
    predicted_loc(i,:) = [xq(image_j(i)),yq(image_i(i))];
    err1(i)=norm(labels_test(i,:)-predicted_loc(i,:));
end
cdfplot(err1);
disp(median(err1)), hold off;
%%
save(['/media/user1/easystore/datasets/aditya_sim/', 'jacobs_spotfi_errors_w_cart.mat'], ...
            'err1', 'act_fin_error')
%% Plot the train and test set and the errors across the test set scattered
figure(4), scatter(labels_train(:,1),labels_train(:,2),'r','filled'), hold on,...
    scatter(labels_test(:,1),labels_test(:,2),'g+'),legend('Training data','testing data'), 
for i=1:4
    scatter(ap{i}(:,1),ap{i}(:,2),'k','filled');
end
hold off;
error_increase = err1 - act_fin_error;
error_locations_cnn = zeros(m,n);
error_locations_spotfi = zeros(m,n);
error_locations_spotfi_cnn_increase = zeros(m,n);
for i=1:n_pts
    error_locations_cnn(li(i),lj(i)) = err1(i);
    error_locations_spotfi(li(i),lj(i)) = act_fin_error(i);
    error_locations_spotfi_cnn_increase(li(i),lj(i)) = error_increase(i);
end
[D1,D2] = meshgrid(d1,d2);
figure(5),subplot(1,2,1), imagesc(d1,d2,error_locations_cnn), axis xy, title('e115 CNN error location map'), colorbar, hold on,
for i=1:4
    scatter(ap{i}(:,1),ap{i}(:,2),'k','filled');
end
hold off;
subplot(1,2,2), plot(labels_fov(:,1),labels_fov(:,2),'r'), plot(labels_nonfov(:,1),labels_nonfov(:,2),'r'),xlim([-9 27]),ylim([-4 12]), title('Test labels Pose Graph'), hold on,
for i=1:4
    scatter(ap{i}(:,1),ap{i}(:,2),'k','filled');
end
hold off;
%%
figure(9),
for i=1:100:length(high_difference_indices)    
    subplot(2,2,1), imagesc(d1,d2,squeeze(geomean(features_w_2d(high_difference_indices(i),:,:,:),2))), axis xy, hold on,title('With offset geomean'),...
        scatter(predicted_loc(high_difference_indices(i),1),predicted_loc(high_difference_indices(i),2),'r*'),...
        scatter(labels_test(high_difference_indices(i),1),labels_test(high_difference_indices(i),2),'g*'),
    for j=1:4
        scatter(ap{j}(:,1),ap{j}(:,2),'k*');
    end
    hold off;
    subplot(2,2,2), imagesc(d1,d2,squeeze(geomean(features_wo_2d(high_difference_indices(i),:,:,:),2))), axis xy, hold on,title('WO offset geomean'),...
        scatter(predicted_loc(high_difference_indices(i),1),predicted_loc(high_difference_indices(i),2),'r*'),...
        scatter(labels_test(high_difference_indices(i),1),labels_test(high_difference_indices(i),2),'g*'),
    for j=1:4
        scatter(ap{j}(:,1),ap{j}(:,2),'k*');
    end
    hold off;
    
    map_sel = outputs(high_difference_indices(i));
    map_xel = map_sel{1,1};
    
    subplot(2,2,3), imagesc(d1,d2,squeeze(geomean(wo_map_xel(:,:,:),2))),axis xy,title('CNN output Offset compensated images geomean'),hold on,...
        scatter(predicted_loc(high_difference_indices(i),1),predicted_loc(high_difference_indices(i),2),'r*'),...
        scatter(labels_test(high_difference_indices(i),1),labels_test(high_difference_indices(i),2),'g*'), hold off;
    subplot(2,2,4), imagesc(d1,d2,squeeze(map_xel(1,:,:))),axis xy,title('CNN prediction'),hold on,...
        scatter(predicted_loc(high_difference_indices(i),1),predicted_loc(high_difference_indices(i),2),'r*'),...
        scatter(labels_test(high_difference_indices(i),1),labels_test(high_difference_indices(i),2),'g*'), hold off;
    %waitforbuttonpress;
end