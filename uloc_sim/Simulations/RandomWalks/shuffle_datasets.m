clear
n_max=25000;
for dataset_idx=2:7
    load(sprintf('Datasets/dataset_d1_d2_%d.mat',dataset_idx));
    load(sprintf('Datasets/ordered_idx_%d.mat',dataset_idx));
    features=features(ordered_idx(1:n_max),:,:,:);
    labels_discrete=labels_discrete(ordered_idx(1:n_max),:);
    labels_gaussian_2d=labels_gaussian_2d(ordered_idx(1:n_max),:,:);
    save(sprintf('Datasets/dataset_d1_d2_%d_shuffled.mat',dataset_idx),'features','labels_discrete','ap','cur_model','labels_gaussian_2d','-v7.3');
end