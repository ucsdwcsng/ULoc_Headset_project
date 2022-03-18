function aoa = get_aoa_of_channels_abs(channels,ap,theta_vals,d_vals,opt)

n_ap=length(ap);
% n_ant=length(ap{1});

% channels_rel = zeros(n_lambda,n_ap,n_ant,n_ap-1);
% features = zeros(n_ap,length(d2),length(d1));
% feature_idx=1;
aoa = zeros(n_ap,1);

for j=1:n_ap
    P = compute_multipath_profile2d_fast_edit(squeeze(channels(:,:,j)),theta_vals,d_vals,opt);
    aoa(j) = get_aoa_for_least_tof(P,d_vals,theta_vals);
end