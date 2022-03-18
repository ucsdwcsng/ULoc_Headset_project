function [features] = generate_features_least_tof_abs(cur_channels,ap,theta_vals,d_vals,d1,d2,opt)

n_ap=size(cur_channels,3);
features = zeros(n_ap*(n_ap-1),length(d2),length(d1));
ap_orient = [1,2,1,2];
for j=1:n_ap
        Prel = compute_multipath_profile2d_fast_edit(squeeze(cur_channels(:,:,j)),theta_vals,d_vals,opt);
        [P,strength] = get_aoa_least_tof1(Prel,d_vals,theta_vals);
        
        ant_pos = abs(ap{j}(1,ap_orient(j))-ap{j}(:,ap_orient(j)));
        channels_new = compute_channels_from_multipath_profile(P,ant_pos,opt.lambda,theta_vals);
        
        Pnew = compute_multipath_profile(channels_new,ant_pos,opt.lambda,theta_vals);
        P_out = strength*convert_multipath_to_2d(Pnew,theta_vals,d1,d2,ap{j});
        P_out = P_out./max(P_out(:));
        
        features(j,:,:) = db(abs(P_out));
end
end
%%
function [aoa,strength] = get_aoa_least_tof1(channel,d_vals,theta_vals)
aoa = zeros(length(theta_vals),1);
channels = abs(channel);
[m,n] = size(channels);
if(m==length(d_vals))
    channels = channels.';
elseif(n~=length(d_vals))
    error('Check your channel matrix dimensions')
end
bw = imregionalmax(db(channels));
[ao,tof] = find(bw);
tof(tof==1) = [];
tof(tof==length(d_vals)) = [];
ao(ao==1) = [];
ao(ao==length(theta_vals)) = [];
if(isempty(tof) || isempty(ao))
     [~,idx] =  max(abs(channel(:)));
     [ao,tof] = ind2sub(size(channels),idx);
end
[~,ind] = min(tof);
aoa(ao(ind)) = 1;
strength = channels(ao(ind),tof(ind));
end    