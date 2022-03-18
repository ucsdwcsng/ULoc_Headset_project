clear
%% Description: 
% The goal of this code is to get channels given a physical space. 
% The definition of physical space can include finite reflectors and 
% bounded obstacles.
% Reflectors are obstacles by default
% Color code: Red is for walls, blue is for obstacles and green is for
% reflectors
% For paths, red path is a non-blocked path and black path is a blocked
% path
%% Define the space
max_x = 20;
max_y = 15;
n_ap = 4;
d1 = 0:0.5:max_x;
d2=0:0.5:max_y;
n_total = 100000;
n_total_per_ap = 20000;
n_ant_per_ap = 3;
features = zeros(n_total,n_ap,length(d2),length(d1));
labels = zeros(n_total,2);

for idx=1:5
    %model=load_model(idx);
    model=load_model(1);
    %% Define the setup
    max_x = max(model.walls(:,1));
    max_y = max(model.walls(:,2));
    ant_sep = min(model.lambda)/2;
    ap{1}=[10,0;
        10+ant_sep,0;
        10+2*ant_sep,0];

    ap{2}=[10,15;
        10+ant_sep,15;
        10+2*ant_sep,15];

    ap{3}=[0,7.5-ant_sep;
            0, 7.5;
            0,7.5+ant_sep;];

    ap{4}=[20,7.5-ant_sep;
            20, 7.5;
            20,7.5+ant_sep;];
    ap_orient=[1,1,2,2];
    
    theta_vals = [-pi/2:0.01:pi/2];
    
    %% Generate the points. Ap 1 antenna 1 serves as the reference
    %features=zeros(n_total,length(model.lambda)*length(ap)*n_ant_per_ap*2);
    
    channels = zeros(length(model.lambda),length(ap),n_ant_per_ap);
    rng;

    start_time = now;
    
    for i=(idx-1)*n_total_per_ap+1:idx*n_total_per_ap
        point_found=false;
        while(~point_found) % Don't allow points inside obstacles
            pos = [rand()*max(model.walls(:,1)),rand()*max(model.walls(:,2))];
            point_found=true;
            for obst_idx = 1:length(model.obstacles)
                if(inpolygon(pos(1),pos(2),model.obstacles{obst_idx}(:,1),model.obstacles{obst_idx}(:,2)))
                    point_found=false;
                end
            end
        end
        for j=1:length(ap)
            for k=1:n_ant_per_ap
                channels(:,j,k)=get_channels_from_model(model,pos,ap{j}(k,:),false);
            end
            channels(:,j,:) = awgn(squeeze(channels(:,j,:)),30);
            P = compute_multipath_profile(squeeze(channels(:,j,:)),ap{j}(:,ap_orient(j)),model.lambda,theta_vals);
            P_out = convert_multipath_to_2d(P,theta_vals,d1,d2,ap{j});
            %channels = channels./repmat(channels(:,1,1),1,length(ap),n_ant_per_ap);
        
            features(i,j,:,:) = abs(P_out);
        end
        labels(i,:) = pos;
        if(mod(i,100)==0)
            disp([i,(now -start_time)*24*60]);
        end
    end
end
%delete(gcp('nocreate'));
%% Convert the labels into one hot encoding
labels_discrete = ceil(labels./2);
n_xlabels = length(unique(labels_discrete(:,1)));
n_ylabels =length(unique(labels_discrete(:,2)));
labels_idx=sub2ind([n_xlabels,n_ylabels],labels_discrete(:,1),labels_discrete(:,2));
labels_one_hot = zeros(size(labels,1),n_xlabels*n_ylabels);
for i=1:size(labels_one_hot,1)
    labels_one_hot(i,labels_idx(i))=1;
    
end
%% Split the dataset into train and test
train_features = features(1:100000,:);
test_features = features(100001:120000,:);
train_labels = labels_one_hot(1:100000,:);
test_labels = labels_one_hot(100001:120000,:);
clear model
for idx=1:6
    model{idx} = load_model(idx);
end
save('dataset_multipath_profiles_one_hot_2m_1model.mat','model','features','labels_one_hot','ap','-v7.3');