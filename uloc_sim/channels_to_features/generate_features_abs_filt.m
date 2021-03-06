function features = generate_features_abs_filt(channels,ap,theta_vals,d_vals,d1,d2,opt)

n_ap=length(ap);
n_ant=size(ap{1},1);

% channels_rel = zeros(n_lambda,n_ap,n_ant,n_ap-1);
features = zeros(n_ap,length(d2),length(d1));
% feature_idx=1;

for j=1:n_ap
    P = compute_multipath_profile2d_fast_edit(squeeze(channels(:,:,j)),theta_vals,d_vals,opt);
    img = abs(P);
    img = (img-min(img(:)))./(max(img(:))-min(img(:)));  
    [thresh, em] = graythresh(img);
    P_new = P.*double(imbinarize(img,thresh*em));
    
    P_out = convert_spotfi_to_2d(P_new,theta_vals,d_vals,d1,d2,ap{j});
    features(j,:,:) = abs(P_out);
end
end