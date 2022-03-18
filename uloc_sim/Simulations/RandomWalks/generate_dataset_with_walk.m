clear
num_pts = 5000;
rel_margin_x = 1/3; % relative size of margin
rel_margin_y = 1/3;
frac = 50;
for dataset_idx=[7]
    storage_path = '/home/yguddeti/aarun/data/p2slam_sim/w_margin/';
%     fn=sprintf('%schannels_revloc_w_margin_%d',storage_path, dataset_idx);
%     load(sprintf('%s.mat',fn));
    fn=sprintf('%slabels_%d',storage_path, dataset_idx);
    
    cur_model=load_model_edit(dataset_idx);
    % Define the setup
    min_x =min(cur_model.walls(:,1));
    min_y =min(cur_model.walls(:,2));
    max_x =max(cur_model.walls(:,1));
    max_y =max(cur_model.walls(:,2));
    margin_x = (max_x - min_x)*rel_margin_x;
    margin_y = (max_y - min_y)*rel_margin_y;
    opt.step_size=0.05;
    opt.walk_pts=100;
    opt.margins=0;
    labels = [rand([num_pts, 1])*(max_x-min_x)*(1-2*rel_margin_x) + min_x + margin_x, ...
                rand([num_pts, 1])*(max_y-min_y)*(1-2*rel_margin_y) + min_y + margin_y];
%     labels_avg = mean(labels, 2);
%     labels_avg_noise = mean(labels_noise, 2);
    walk_pts=zeros(size(labels,1)/frac,opt.walk_pts);
%     walk_pts_noisy=zeros(size(labels,1),opt.walk_pts);
    for i=1:size(labels,1)/frac
        if mod(i, 10) == 0
            disp(i)
        end
        walk_idx=get_path_from_dataset(labels , 'none', opt);
%         walk_idx_noisy=get_path_from_dataset(labels_avg_noise, 'none', opt);
        walk_pts(i,:)=walk_idx';
%         walk_pts_noisy(i,:)=walk_idx_noisy';
    end
%     save(sprintf('%s_with_walk_100.mat',fn),'walk_pts', 'walk_pts_noisy', 'ap', ...
%                     'cur_model', 'channels_w_offset', 'channels_wo_offset', ...
%                     'channels_w_offset_noisy', 'channels_wo_offset_noisy', ...
%                     'labels', 'labels_noise', '-v6');
    save(sprintf('%s_with_walk_100.mat',fn),'walk_pts', 'cur_model', 'labels', '-v6');
end