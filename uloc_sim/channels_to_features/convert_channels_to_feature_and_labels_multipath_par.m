clear;

for j=9
    disp(j);
    cur_model = load_model_edit(j);
    ant_sep = min(cur_model.lambda)/2;
    ap{1}=[10,0;
        10+ant_sep,0;
        10+2*ant_sep,0];

    ap{2} = [10,15;
        10+ant_sep,15;
        10+ant_sep*2,15];

    ap{3} = [0,7.5-ant_sep;
            0,7.5;
            0,7.5+ant_sep];

    ap{4} = [20,7.5-ant_sep;
            20,7.5;
            20,7.5-2*ant_sep]; 
    input_grid_size=0.5;  %Grid size for input Gaussian
    output_grid_size=0.5;   %Grid size for output Gaussian
    output_sigma = 1; %Sigma for Gaussian output image
    ap_orient=[1,1,2,2];% 1 if ap is along x, 2 otherwise
    theta_vals = [-pi/2:0.01:pi/2]; %Theta values to consider for multipath profiles
    d_vals = -15:0.1:15; %D vals to consider for distance profiles
    n_points=40000;
    n_lambda=length(cur_model.lambda);
    n_ap=length(ap);
    n_ant=length(ap{1});

    n_features=n_ap*n_ap; % One distance profile per antenna, one angle profile per ap
    max_x =max(cur_model.walls(:,1));
    max_y =max(cur_model.walls(:,2));
    d1 = 0:input_grid_size:max_x;
    d2 = 0:input_grid_size:max_y;
    
%     labels_gaussian = zeros(size(labels,1),n_xlabels*n_ylabels);
%     map_X = repmat((1:n_xlabels)'*output_grid_size,1,n_ylabels);
%     map_Y = repmat((1:n_ylabels)*output_grid_size,n_xlabels,1);

    channels_rel= zeros(length(cur_model.lambda),n_ap,n_ant,n_ap-1);
    opt.freq = 3e8./(cur_model.lambda);
    opt.ant_sep = abs(ap{1}(2,1)-ap{1}(1,1));
    
    features = zeros(n_points,n_ap*n_ap,length(d2),length(d1));
    labels_discrete = zeros(n_points,2);
    i=0;
    while( i<=n_points)        
        point_found = false;
        while(~point_found) % Don't allow points inside obstacles
            pos = [rand()*max(cur_model.walls(:,1)),rand()*max(cur_model.walls(:,2))];
            point_found=true;
            for obst_idx = 1:length(cur_model.obstacles)
                if(inpolygon(pos(1),pos(2),cur_model.obstacles{obst_idx}(:,1),cur_model.obstacles{obst_idx}(:,2)))
                    point_found=false;
                end
            end
        end
        i=i+1;
        % Get features
        features(i,:,:,:) = generate_features_rel_only(pos,ap,theta_vals,d_vals,d1,d2,opt,cur_model);
        % Get labels
        labels_discrete(i,:) = pos;
    %         d = (map_X-labels(i,1)).^2+(map_Y-labels(i,2)).^2;
    %         cur_gaussian = exp(-d/output_sigma/output_sigma)*1/sqrt(2*pi)/output_sigma;
    %         labels_gaussian(i,:)=cur_gaussian(:);
        if(mod(i,1000)==0)
            disp(i);
        end
        
    end
    label_std = 1/sqrt(2*pi);%std(labels_gaussian_all(:));
    feature_mean =zeros(1,size(features,2));
    feature_std = zeros(1,size(features,2));
    for k=1:length(feature_mean)
        dset = features(:,k,:,:);
        feature_mean(k) = mean(dset(:));
        feature_std(k) = std(dset(:));
    end
    features = (features-repmat(feature_mean,size(features,1),1,size(features,3),size(features,4)))./repmat(feature_std,size(features,1),1,size(features,3),size(features,4));
    labels_gaussian_2d=get_gaussian_labels1(labels_discrete,output_grid_size,2);
    stri = ['datasets/dataset_d1_d2_',num2str(j),'.mat'];
    save(stri,'features','labels_discrete','labels_gaussian_2d','cur_model','ap','-v7.3');
    clear
end



