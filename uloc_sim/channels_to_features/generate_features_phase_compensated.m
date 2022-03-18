function features = generate_features_phase_compensated(channels,ap,d_vals,theta_vals,d1,d2,opt)

n_ap=length(ap);
% n_ant=length(ap{1});

% channels_rel = zeros(n_lambda,n_ap,n_ant,n_ap-1);
features = zeros(n_ap,length(d2),length(d1));
feature = zeros(n_ap,length(d_vals),length(theta_vals));
% feature_idx=1;

for j=1:n_ap
    P = compute_multipath_profile2d_fast_edit(squeeze(channels(:,j,:)),theta_vals,d_vals,opt);
    feature(j,:,:) = abs(P);
    offset =  aoa_tof_intersect(feature,d_vals,theta_vals,ap);
    P = compute_multipath_profile2d_offset_compensation(squeeze(channels(:,j,:)),theta_vals,d_vals,opt,offset);
    P_out = convert_spotfi_to_2d(P,theta_vals,d_vals,d1,d2,ap{j});
    features(j,:,:) = abs(P_out);
end
end
