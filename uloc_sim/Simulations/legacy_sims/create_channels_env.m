clearvars
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
n_ap = 4;
% d1 = 0:0.1:max_x;
% d2=0:0.1:max_y;
% d1 = -2:0.1:2;
% d2 = -2:0.1:2;
n_total_per_ap = 10000;
n_ant_per_ap = 4;
labels = zeros(n_total_per_ap,2);
for idx=2:6
    disp(idx);
    cur_model=load_model_edit(idx);
    %model=load_model(1);
    %% Define the setup
    min_x =min(cur_model.walls(:,1));
    min_y =min(cur_model.walls(:,2));
    max_x =max(cur_model.walls(:,1));
    max_y =max(cur_model.walls(:,2));
    ant_sep = min(cur_model.lambda)/2;
%     dist_ant = [-1.5*ant_sep,-0.5*ant_sep,0.5*ant_sep,1.5*ant_sep];
%     ap{1} = [5*ones(4,1)+dist_ant',zeros(4,1)];
%     ap{2} = [5*ones(4,1)+dist_ant',8*ones(4,1)];
%     ap{3} = [zeros(4,1),4*ones(4,1)+dist_ant'];
%     ap{4} = [10*ones(4,1),4*ones(4,1)+dist_ant'];
%     ap{1} = [0.643,2.5+1.5*ant_sep;...
%         0.643,2.5+0.5*ant_sep;...
%         0.643,2.5-0.5*ant_sep;...
%         0.643,2.5-1.5*ant_sep];
%     ap{2} = [6.65-1.5*ant_sep,0.055;...
%         6.65-0.5*ant_sep,0.055;...
%         6.65+0.5*ant_sep,0.055;...
%         6.65+1.5*ant_sep,0.055];
%     ap{3} = [12.279,2.78-1.5*ant_sep;...
%         12.279,2.78-0.5*ant_sep;...
%         12.279,2.78+0.5*ant_sep;...
%         12.279,2.78+1.5*ant_sep];
%     ap{4} = [4.262+1.5*ant_sep,4.611;...
%         4.262+0.5*ant_sep,4.611;...
%         4.262-0.5*ant_sep,4.611;...
%         4.262-1.5*ant_sep,4.611];
ap{1} = [4.2-1.5*ant_sep,0;...
    4.2-0.5*ant_sep,0;...
    4.2+0.5*ant_sep,0;...
    4.2+1.5*ant_sep,0];
ap{2} = [8.36,2-1.5*ant_sep;...
    8.36,2-0.5*ant_sep;...
    8.36,2+0.5*ant_sep;...
    8.36,2+1.5*ant_sep];
ap{3} = [3+1.5*ant_sep,4.15;...
    3+0.5*ant_sep,4.15;...
    3-0.5*ant_sep,4.15;...
    3-1.5*ant_sep,4.15];
ap{4} = [0,2+1.5*ant_sep;...
    0,2+0.5*ant_sep;...labels
    0,2-0.5*ant_sep;...
    0,2-1.5*ant_sep];
    Ts = cur_model.Ts;
%     ap{1}=[,0;
%         10+ant_sep,0;
%         10+2*ant_sep,0];
% 
%     ap{2}=[10,15;
%         10+ant_sep,15;
%         10+2*ant_sep,15];
% 
%     ap{3}=[0,7.5-ant_sep;
%             0, 7.5;
%             0,7.5+ant_sep;];
% 
%     ap{4}=[20,7.5-ant_sep;
%             20, 7.5;
%             20,7.5+ant_sep;];
    ap_orient=[2,1,2,1];
    
    theta_vals = [-pi/2:0.01:pi/2];
    
    %% Generate the points. Ap 1 antenna 1 serves as the reference
    %features=zeros(n_total,length(model.lambda)*length(ap)*n_ant_per_ap*2);
    
    channels_w_offset = zeros(n_total_per_ap,length(cur_model.lambda),length(ap),n_ant_per_ap);
    channels_wo_offset = zeros(n_total_per_ap,length(cur_model.lambda),length(ap),n_ant_per_ap);
    offsets = zeros(n_total_per_ap,4);
    rng;

    start_time = now;
    
    parfor i=1:n_total_per_ap
%         point_found=false;
%         while(~point_found) % Don't allow points inside obstacles
       pos = [rand()*(max_x-min_x)+min_x,rand()*(max_y-min_y)+min_y];
%             point_found=true;
%             for obst_idx = 1:length(model.obstacles)
%                 if(inpolygon(pos(1),pos(2),model.obstacles{obst_idx}(:,1),model.obstacles{obst_idx}(:,2)))
%                     point_found=false;
%                 end
%             end
%         end
        [channels_wo_offset(i,:,:,:),channels_w_offset(i,:,:,:),offset(i,:)] = call_get_channels_from_model(cur_model,cur_model,pos,ap);
        labels(i,:) = pos;
        if(mod(i,1000)==0)
            disp(i);
        end
    end
    save(sprintf('/media/user1/easystore/datasets/pix2pix/channels/channels_pix2pix_9x5_%d.mat',idx),'cur_model','channels_wo_offset','channels_w_offset','offset','labels','ap','-v7.3');
end
%delete(gcp('nocreate'));
