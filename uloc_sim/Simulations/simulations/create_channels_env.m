clc
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
% max_x = 13;
% max_y = 5;
disp('Simulation starting')
n_ap = 4;
load('noise_stats.mat')
storage_path = '/home/aarun/Research/Data/p2slam_sim/';

% storage_path = './';
% d1 = 0:0.1:max_x;
% d2=0:0.1:max_y;
% d1 = -2:0.1:2;
% d2 = -2:0.1:2;
% n_total_per_ap = 5000;
n_rx_ant = 4;
n_tx_ant = 4;
rel_pos = [-1.5;-0.5;0.5;1.5];

% rel_margin_x = 1/3; % relative size of margin
% rel_margin_y = 1/3;
% flags
add_pos_noise = true;
save_file = true;
%%
for idx=[7]
    fn=sprintf('%slabels_%d_with_walk_100',storage_path, idx);
    load(sprintf('%s.mat',fn));
    labels_single = labels;
%     labels = zeros(n_total_per_ap, n_tx_ant, 2);
    labels_noise = zeros(size(walk_pts, 1),size(walk_pts, 2), ...
                         n_tx_ant, 2);
%     orients = zeros(n_total_per_ap, 1);
    storage_path = '/home/yguddeti/aarun/data/p2slam_sim/w_margin/';

    disp(idx);
    cur_model=load_model_edit(idx);
    % Define the setup
%     min_x =min(cur_model.walls(:,1));
%     min_y =min(cur_model.walls(:,2));
%     max_x =max(cur_model.walls(:,1));
%     max_y =max(cur_model.walls(:,2));
%     margin_x = (max_x - min_x)*rel_margin_x;
%     margin_y = (max_y - min_y)*rel_margin_y;
    ant_sep = min(cur_model.lambda)/2;

    ap{1} = [0,2+1.5*ant_sep;...
            0,2+0.5*ant_sep;...
            0,2-0.5*ant_sep;...
            0,2-1.5*ant_sep];
    ap{2} = [16.8,7.6-1.5*ant_sep;...
            16.8,7.6-0.5*ant_sep;...
            16.8,7.6+0.5*ant_sep;...
            16.8,7.6+1.5*ant_sep];
    ap{3} = [6.4+1.5*ant_sep,7.6;...
            6.4+0.5*ant_sep,7.6;...
            6.4-0.5*ant_sep,7.6;...
            6.4-1.5*ant_sep,7.6];
    ap{4} = [12.4-1.5*ant_sep,0;...
            12.4-0.5*ant_sep,0;...
            12.4+0.5*ant_sep,0;...
            12.4+1.5*ant_sep,0];

    Ts = cur_model.Ts;
    
    theta_vals = [-pi/2:0.01:pi/2];
    
    %Generate the points. Ap 1 antenna 1 serves as the reference
    %features=zeros(n_total,length(model.lambda)*length(ap)*n_rx_ant*2);
    
    channels_w_offset = zeros(size(walk_pts, 1),size(walk_pts, 2), ...
                                length(cur_model.lambda), ...
                                length(ap),n_rx_ant, n_tx_ant);
    channels_wo_offset = zeros(size(walk_pts, 1),size(walk_pts, 2), ...
                                length(cur_model.lambda), ...
                                length(ap),n_rx_ant, n_tx_ant);
    if add_pos_noise
        channels_w_offset_noisy = zeros(size(walk_pts, 1),size(walk_pts, 2), ...
                                    length(cur_model.lambda), ...
                                    length(ap),n_rx_ant, n_tx_ant);
        channels_wo_offset_noisy = zeros(size(walk_pts, 1),size(walk_pts, 2), ...
                                    length(cur_model.lambda), ...
                                    length(ap),n_rx_ant, n_tx_ant);
    end
    offset = zeros(size(walk_pts, 1),size(walk_pts, 2), ...
                    n_rx_ant, n_tx_ant);
    
    for row=1:size(walk_pts, 1)
        cur_wp = walk_pts(row, :);
        heading = labels_single(cur_wp(2:end), :) ...
                        - labels_single(cur_wp(1:end-1), :);
        heading = atan2(heading(:, 2), heading(:, 1));
        heading = [heading; heading(end)];
        for col=1:size(walk_pts, 2)
%             rng(row+idx*n_total_per_ap);
    %         pos = [rand()*(max_x-min_x)*(1-2*rel_margin_x) + min_x + margin_x, ...
    %                 rand()*(max_y-min_y)*(1-2*rel_margin_y) + min_y + margin_y];
    %         th = [rand()*pi];
            pos = labels_single(walk_pts(row, col), :);
            th = heading(col);
            delta_loc = ant_sep*[cos(th), sin(th)];
            pos_4 = pos + delta_loc.*rel_pos; % hardcoded or 4 Tx antennas
            cur_offset = (rand(1,4)-0.25)*6.*cur_model.Ts;
            if (add_pos_noise)
                noise_x = noise_std(1)*rand() + noise_mean(1);
                noise_y = noise_std(2)*rand() + noise_mean(2);
                pos_noise = pos_4 + [noise_x, noise_y];
                for j=1:n_tx_ant
                    [channels_wo_offset_noisy(row, col,:,:,:,j), ...
                    channels_w_offset_noisy(row,col,:,:,:,j), ...
                    offset(row,col,:,j)] = ...
                        call_get_channels_from_model(cur_model,cur_model,pos_noise(j,:),ap, cur_offset);
                end
                labels_noise(row,col, :, :) = pos_noise;
            end
            for j = 1:n_tx_ant
                [channels_wo_offset(row,col,:,:,:, j), ...
                channels_w_offset(row,col,:,:,:, j), ...
                offset(row,col,:, j)] = ...
                            call_get_channels_from_model(cur_model,cur_model,pos_4(j, :),ap, cur_offset);
            end
%             labels(row,:, :) = pos_4;
    %         orients(i, 1) = th;

        end
        if(mod(row,100) == 0)
            disp(row);
        end
    end
    file_name = [storage_path, ...
                sprintf('channels_revloc_w_margin_%d_with_walk_100.mat',idx)];
    if save_file
        if (add_pos_noise)
            save(file_name,'cur_model','channels_wo_offset','channels_w_offset', ...
                         'channels_wo_offset_noisy','channels_w_offset_noisy', ...
                         'offset','labels', 'labels_noise', 'walk_pts', 'ap','-v7.3');
        else
            save(file_name,'cur_model','channels_wo_offset','channels_w_offset', ...
                         'offset','labels', 'ap','-v7.3');
        end
    end
end
