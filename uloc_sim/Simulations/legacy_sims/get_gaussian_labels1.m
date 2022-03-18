function labels_gaussian = get_gaussian_labels1(labels,grid_size,sigma)
    labels_discrete = round(labels./grid_size);
    n_xlabels = length(unique(labels_discrete(:,1)));
    n_ylabels =length(unique(labels_discrete(:,2)));
    labels_gaussian = zeros(size(labels,1),n_ylabels,n_xlabels);
    map_X = repmat((1:n_xlabels)*grid_size,n_ylabels,1);
    map_Y = repmat((1:n_ylabels)'*grid_size,1,n_xlabels);
    n_points = size(labels,1);
        for i=1:n_points
            d = (map_X-labels(i,1)).^2+(map_Y-labels(i,2)).^2;
            cur_gaussian = exp(-d/sigma/sigma);%*1/sqrt(2*pi)/output_sigma;        
            labels_gaussian(i,:)=cur_gaussian(:);
        end
end



