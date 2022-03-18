function[reg_max] = get_aoa_tof_features(channels,ap_pos,theta_vals,d_vals,opt)

% features = generate_aoa_tof_features(channels,ap_pos,theta_vals,d_vals,opt);
n_ap=length(ap_pos);
% n_ant=length(ap{1});

% channels_rel = zeros(n_lambda,n_ap,n_ant,n_ap-1);
features = zeros(n_ap,length(theta_vals),length(d_vals));
% feature_idx=1;

for j=1:n_ap
    P = compute_multipath_profile2d_fast_edit(squeeze(channels(:,j,:)),theta_vals,d_vals,opt);
    features(j,:,:) = abs(P);
end

reg_max = zeros(size(features,1),2);
for i=1:size(features,1)
    reg_max(i,:) = aoa_tof_max(squeeze(features(i,:,:)),d_vals,theta_vals);
end
end