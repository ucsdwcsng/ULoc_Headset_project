%clear
%load dataset_multipath_profiles.mat
% This script will convert the dataset into one hot encoding
grid_size = 2;
labels_discrete = ceil(labels./grid_size);
n_xlabels = length(unique(labels_discrete(:,1)));
n_ylabels =length(unique(labels_discrete(:,2)));
sigma=1;
labels_gaussian_2d=zeros(size(labels,1),n_xlabels,n_ylabels);

labels_gaussian = zeros(size(labels,1),n_xlabels*n_ylabels);
for i=1:size(labels_discrete,1)
   map_X = repmat((1:n_xlabels)'*grid_size,1,n_ylabels);
   map_Y = repmat((1:n_ylabels)*grid_size,n_xlabels,1);
   d = (map_X-labels(i,1)).^2+(map_Y-labels(i,2)).^2;
   cur_gaussian = exp(-d/sigma/sigma)*1/sqrt(2*pi)/sigma;
   labels_gaussian(i,:)=cur_gaussian(:);
   % curr_gaussian = 
    
end
    