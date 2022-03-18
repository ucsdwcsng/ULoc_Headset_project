clear
% channel
for dataset_idx=7
    load(sprintf('Datasets/dataset_d1_d2_%d.mat',dataset_idx));
    ordered_idx=zeros(size(labels_discrete,1),1);
    selected_idx=false(size(labels_discrete,1),1);
    dist_sep=0.05;
    ordered_idx(1)=1;
    selected_idx(1)=true;
    for i=2:size(ordered_idx)
        cur_pt=labels_discrete(ordered_idx(i-1),:);
        available_idx=find((~selected_idx) & abs(labels_discrete(:,1)-cur_pt(1))<dist_sep & abs(labels_discrete(:,2)-cur_pt(2))<dist_sep);
        if(isempty(available_idx))
            disp(i);
            break;
        end
        pick_idx = randi(length(available_idx));
        selected_idx(available_idx(pick_idx))=true;
        ordered_idx(i)=available_idx(pick_idx);
        %disp(labels_discrete(ordered_idx(i),:));
    end    
    ordered_idx=ordered_idx(1:i-1);
    save(sprintf('Datasets/ordered_idx_%d',dataset_idx),'ordered_idx');
end