function features = generate_features_rel_only(pos,ap,theta_vals,d_vals,d1,d2,opt,model)
n_lambda=size(model.lambda,2);
n_ap=4;
n_ant=length(ap{1});
channels = zeros(n_lambda,n_ap,n_ant);
channels_rel = zeros(n_lambda,n_ap,n_ant,n_ap-1);

features = zeros(n_ap*n_ap,length(d2),length(d1));
for j=1:n_ap
    for k=1:n_ant
        channels(:,j,k)=get_channels_from_model_edit(model,pos,ap{j}(k,:),false);
        channels(:,j,k)=awgn(squeeze(channels(:,j,k)),30,0,'db');
    end
end
feature_idx=1;
for j=1:n_ap
    for k=1:n_ant
        ind = 1:1:n_ap;
        ind(j) = [];
        for l = 1:n_ap-1
            channels_rel(:,j,k,l) = channels(:,j,k).*conj(channels(:,ind(l),1));
        end
    end

    for l = 1:n_ap-1

        Prel =  compute_multipath_profile2d_fast_edit(squeeze(channels_rel(:,j,:,l)),theta_vals,d_vals,opt);
        Prel_out = convert_relative_spotfi_to_2d_edit(Prel,ap{j},ap{ind(l)},theta_vals,d_vals,d1,d2);
        features(feature_idx,:,:) = db(abs(Prel_out));
        feature_idx = feature_idx+1;

    end
    P = compute_multipath_profile2d_fast_edit(squeeze(channels(:,j,:)),theta_vals,d_vals,opt);
    P_out = convert_spotfi_to_2d(P,theta_vals,d_vals,d1,d2,ap{j});

    features(feature_idx,:,:) = db(abs(P_out));
    
    feature_idx = feature_idx+1;
end
    