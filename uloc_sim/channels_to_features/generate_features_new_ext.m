function features = generate_features_new_ext(channels,ap,theta_vals,d_vals,d1,d2,opt,model)

n_lambda=size(model.lambda,2);
n_ap=4;
n_ant=length(ap{1});

channels_rel = zeros(n_lambda,n_ap,n_ant,(n_ap-1)*n_ant);
features = zeros(n_ap*(n_ap-1)*n_ant,length(d2),length(d1));
feature_idx=1;

for j=1:n_ap
    for k=1:n_ant
        for ref_ant = 1:n_ant
            ind = 1:1:n_ap;
            ind(j) = [];
            for l = 1:n_ap-1
                channels_rel(:,j,k,(ref_ant-1)*n_ant+l) = channels(:,j,k).*conj(channels(:,ind(l),ref_ant));
            end
        end
    end
    for ref_ant = 1:n_ant
        for l = 1:n_ap-1
            Prel =  compute_multipath_profile2d_fast_edit(squeeze(channels_rel(:,j,:,(ref_ant-1)*n_ant+l)),theta_vals,d_vals,opt);
            Prel_out = convert_relative_spotfi_to_2d_edit(Prel,ap{j},ap{ind(l)},theta_vals,d_vals,d1,d2);
            features(feature_idx,:,:) = db(abs(Prel_out));
            feature_idx = feature_idx+1;

        end
    end
end