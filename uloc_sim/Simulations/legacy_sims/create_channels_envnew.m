% NOTE: changed AWGN function because of some matlab license issue
clear;
for j=1:6
    disp(j);
    cur_model = load_model_edit(j);
    ant_sep = min(cur_model.lambda)/2;
%     dist_ant = (-1.5:1:1.5)'*ant_sep;
%     ap{1}(:,1) = 5+dist_ant;
%     ap{1}(:,2) = 0;
% 
%     ap{2}(:,1) = 5+dist_ant;
%     ap{2}(:,2) = 8;
    
% (-1.357,1.5)
% (4.65,-0.945)
% (10.279,1.78)
% (2.262,3.611)

%     ap{1} = [2-1.357,3.5+1.5*ant_sep;...
%         2-1.357,3.5+0.5*ant_sep;...
%         2-1.357,3.5-0.5*ant_sep;...
%         2-1.357,3.5-1.5*ant_sep];
%     ap{2} = [6.65-1.5*ant_sep,2-0.945;...
%         6.65-0.5*ant_sep,2-0.945;...
%         6.65+0.5*ant_sep,2-0.945;...
%         6.65+1.5*ant_sep,2-0.945];
%     ap{3} = [12.279,3.78-1.5*ant_sep;...
%         12.279,3.78-0.5*ant_sep;...
%         12.279,3.78+0.5*ant_sep;...
%         12.279,3.78+1.5*ant_sep];
%     ap{4} = [4.262+1.5*ant_sep,5.611;...
%         4.262+0.5*ant_sep,5.611;...
%         4.262-0.5*ant_sep,5.611;...
%         4.262-1.5*ant_sep,5.611];

%     ap{4} = [-0.012+1.5*ant_sep,0.691;...
%         -0.012+0.5*ant_sep,0.691;...
%         -0.012-0.5*ant_sep,0.691;...
%         -0.012-1.5*ant_sep,0.691];
% 
%     ap{3} = [8.918+1.5*ant_sep,5.455;...
%         8.918+0.5*ant_sep,5.455;...
%         8.918-0.5*ant_sep,5.455;...
%         8.918-1.5*ant_sep,5.455];
% 
%     ap{2} = [8.830-1.5*ant_sep,0.686;...
%         8.830-0.5*ant_sep,0.686;...
%         8.830+0.5*ant_sep,0.686;...
%         8.830+1.5*ant_sep,0.686];
% 
%     ap{1} = [0.165-1.5*ant_sep,5.527;...
%         0.165-0.5*ant_sep,5.527;...
%         0.165+0.5*ant_sep,5.527;...
%         0.165+1.5*ant_sep,5.527];
    %clear;
    
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
    0,2+0.5*ant_sep;...
    0,2-0.5*ant_sep;...
    0,2-1.5*ant_sep];

% comment: Aditya: why redefine these here? 
%     freq = double(5e9 + 5*153*1e6) + [-117:-1,1:117].*80.*1e6./256;
%     lambda = 3e8./freq;
%     disp(j);
%     cur_model = load_model_edit2(j);
%     cur_model.lambda = lambda;
%     ant_sep = min(cur_model.lambda)/2;
%     dist_ant = (-4:5)'*ant_sep;
%     ap = ap_pos_corrected;
    min_x =min(cur_model.walls(:,1));
    min_y =min(cur_model.walls(:,2));
    max_x =max(cur_model.walls(:,1));
    max_y =max(cur_model.walls(:,2));

    n_points = 10000;
%     labels_gaussian = zeros(size(labels,1),n_xlabels*n_ylabels);
%     map_X = repmat((1:n_xlabels)'*output_grid_size,1,n_ylabels);
%     map_Y = repmat((1:n_ylabels)*output_grid_size,n_xlabels,1);
    n_lambda=length(cur_model.lambda);
    n_ap=length(ap);
    n_ant=length(ap{1});
    channel = zeros(n_points,n_lambda,n_ap,n_ant);

    opt.freq = 3e8./(cur_model.lambda);
    opt.ant_sep = abs(ap{1}(2,1)-ap{1}(1,1));
    label_discrete = zeros(n_points,2);
    i=0;
    offsets = zeros(n_points,n_ap);
    for i=1:n_points    

        pos = [rand()*(max_x-min_x)+min_x,rand()*(max_y-min_y)+min_y];
        for l=1:n_ap
            offset = rand();
            parfor k=1:n_ant
                channels(i,:,l,k)=squeeze(get_channels_from_model_edit(cur_model,pos,ap{l}(k,:),false,offset));
                channels(i,:,l,k)=awgn(squeeze(channels(i,:,l,k)),20);
            end
            offsets(i,l) = offset;
        end
        
        label_discrete(i,:) = pos;
    %         d = (map_X-labels(i,1)).^2+(map_Y-labels(i,2)).^2;
    %         cur_gaussian = exp(-d/output_sigma/output_sigma)*1/sqrt(2*pi)/output_sigma;
    %         labels_gaussian(i,:)=cur_gaussian(:);
        if(mod(i,1000)==0)
            disp(i);
        end 
    end   
    channels(j*n_points + [1:n_points],:,:,:) = channel;
    labels_discrete(j*n_points + [1:n_points],:) = label_discrete;
end
% stri = ['datasets/channel20_10_',num2str(j),'.mat'];
save('/home/aarun/Research/Data/deep_loc/data/sim_data/data_1.mat','channels','labels_discrete','cur_model','ap','-v7.3');
clear
%% display the models:
for i = 1:6
    display_model(load_model_edit1(i))
    waitforbuttonpress
end

