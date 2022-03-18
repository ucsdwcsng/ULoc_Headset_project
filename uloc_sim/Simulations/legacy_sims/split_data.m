clear
fn='dataset_distance_aoa_2';
load(sprintf('%s.mat',fn));
n_split=10;
points_per_split = size(features,1)/n_split;
features_all = features;
labels_gaussian_all = labels_gaussian;
labels_discrete = ceil(labels./2);
n_xlabels = length(unique(labels_discrete(:,1)));
n_ylabels =length(unique(labels_discrete(:,2)));
labels_idx=sub2ind([n_xlabels,n_ylabels],labels_discrete(:,1),labels_discrete(:,2));
labels_one_hot = zeros(size(labels,1),n_xlabels*n_ylabels);
labels_gaussian_2d_all=get_gaussian_labels(labels,0.5,2);
for i=1:size(labels_one_hot,1)
    labels_one_hot(i,labels_idx(i))=1;
    
end
clear features labels_gaussian labels labels_idx
for i=1:n_split
    start_idx = floor((i-1)*points_per_split)+1;
    end_idx = floor(i*points_per_split);
    features = features_all(start_idx:end_idx,:,:,:);
    if(i==1)
        label_std = 1/sqrt(2*pi);%std(labels_gaussian_all(:));
        feature_mean =zeros(1,size(features,2));
        feature_std = zeros(1,size(features,2));
        for j=1:length(feature_mean)
            dset = features(:,j,:,:);
            feature_mean(j) = mean(dset(:));
            feature_std(j) = std(dset(:));
        end
            
    end
    features = (features-repmat(feature_mean,size(features,1),1,size(features,3),size(features,4)))./repmat(feature_std,size(features,1),1,size(features,3),size(features,4));
    labels_gaussian = labels_gaussian_all(start_idx:end_idx,:)/label_std;
    labels_gaussian_2d=labels_gaussian_2d_all(start_idx:end_idx,:,:);
    labels_sum=squeeze(mean(features,2));
    save(sprintf('%s_2dgauss_%d.mat',fn,i),'features','labels_gaussian_2d','labels_one_hot','-v7.3');
end