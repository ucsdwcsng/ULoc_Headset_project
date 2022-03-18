function features = generate_features_least_tof_rel1(cur_channels,ap,theta_vals,d_vals,d1,d2,opt)
n_lambda=size(cur_channels,1);
n_ap=size(cur_channels,2);
n_ant=size(cur_channels,3);
channels = zeros(n_lambda,n_ap,n_ant,n_ap-1);
channels_rel = zeros(n_lambda,n_ap,n_ant,n_ap-1);
% channel_rel = zeros(n_lambda,n_ap,n_ant,n_ap-1);
features = zeros(n_ap*(n_ap-1),length(d2),length(d1));

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
        Prel =  compute_spotfi_profile_edit2(squeeze(channels_rel(:,j,:,l)),theta_vals,d_vals,opt);
%         Prel1 = compute_multipath_profile2d_fast_edit(squeeze(channels_rel(:,j,:,l)),theta_vals,d_vals,opt);
        P = get_aoa_least_tof1(Prel,d_vals);
        P_out = convert_multipath_to_2d(P,theta_vals,d1,d2,ap{j}(:,:));
        features(feature_idx,:,:) = db(abs(P_out));
        feature_idx = feature_idx + 1;
    end
end
end

% function [aoa] = get_aoa_least_tof(channel,d_vals)
% channels = abs(channel);
% [m,n] = size(channels);
% if(m==length(d_vals))
%     channels = channels.';
% elseif(n~=length(d_vals))
%     error('Check your channel matrix dimensions')
% end
% aoa_profile = mean(channels,2);
% [~,inds] = findpeaks(aoa_profile);
% if(inds)
%     aoa_wanted = channels(inds,:);
%     [~,tof_inds] = max(aoa_wanted,[],2);
%     [tof_ind,~]=min(tof_inds,[],1);
%     aoa = channel(:,tof_ind(1));
% 
% else
%     % Peaks are at the edges have been neglected by the findpeaks algorithm
%     inds = [1,length(aoa_profile)]; % assign peak indices as required indices
%     aoa_wanted = channels(inds,:);
%     [~,tof_inds] = max(aoa_wanted,[],2);
%     [tof_ind,~]=min(tof_inds,[],1);
%     aoa = channel(:,tof_ind(1));
% end
% end    
function [aoa] = get_aoa_least_tof1(channel,d_vals)
channels = abs(channel);
[m,n] = size(channels);
if(m==length(d_vals))
    channels = channels.';
elseif(n~=length(d_vals))
    error('Check your channel matrix dimensions')
end
bw = imregionalmax(db(channels));
[test,tof] = find(bw);
tof(tof==1) = [];
tof(tof==length(d_vals)) = [];
if(isempty(tof))
     [~,idx] =  max(abs(channel(:)));
     [test,tof] = ind2sub(size(channels),idx);
end
aoa = channel(:,min(tof));
end   