clear
close all;
%% Load Data
% load('channels_Feb23_calibrated.mat')
% load('ap_cli_pos_real_23Feb.mat')

for j=1:7
    ant_sep = 0.026;
    ap{4} = [-0.012+1.5*ant_sep,0.691;...
        -0.012+0.5*ant_sep,0.691;...
        -0.012-0.5*ant_sep,0.691;...
        -0.012-1.5*ant_sep,0.691];

    ap{3} = [8.918+1.5*ant_sep,5.455;...
        8.918+0.5*ant_sep,5.455;...
        8.918-0.5*ant_sep,5.455;...
        8.918-1.5*ant_sep,5.455];

    ap{2} = [8.830-1.5*ant_sep,0.686;...
        8.830-0.5*ant_sep,0.686;...
        8.830+0.5*ant_sep,0.686;...
        8.830+1.5*ant_sep,0.686];

    ap{1} = [0.165-1.5*ant_sep,5.527;...
        0.165-0.5*ant_sep,5.527;...
        0.165+0.5*ant_sep,5.527;...
        0.165+1.5*ant_sep,5.527];
    %clear;
    freq = double(5e9 + 5*153*1e6) + [-117:-1,1:117].*80.*1e6./256;
    lambda = 3e8./freq;
    disp(j);
    cur_model = load_model_edit(j);
    cur_model.lambda = lambda;
%     ant_sep = min(cur_model.lambda)/2;
%     dist_ant = (-4:5)'*ant_sep;
%     ap = ap_pos_corrected;

    max_x =max(cur_model.walls(:,1));
    max_y =max(cur_model.walls(:,2));

    n_points = 10000;
%     labels_gaussian = zeros(size(labels,1),n_xlabels*n_ylabels);
%     map_X = repmat((1:n_xlabels)'*output_grid_size,1,n_ylabels);
%     map_Y = repmat((1:n_ylabels)*output_grid_size,n_xlabels,1);
    n_lambda=length(lambda);
    n_ap=length(ap);
    n_ant=length(ap{1});
    channels = zeros(n_points,n_lambda,n_ap,n_ant);

    opt.freq = 3e8./(cur_model.lambda);
    opt.ant_sep = abs(ap{1}(2,1)-ap{1}(1,1));
    labels_discrete = zeros(n_points,2);
    i=0;
    offsets = zeros(n_points,n_ap);
    for i=1:n_points    

        pos = [rand()*max_x,rand()*max_y];
        for l=1:n_ap
            offset = rand()*3.75;
            parfor k=1:n_ant
                channels(i,:,l,k)=squeeze(get_channels_from_model_edit(cur_model,pos,ap{l}(k,:),false,offset));
                channels(i,:,l,k)=awgn(squeeze(channels(i,:,l,k)),20,0,'db');
            end
            offsets(i,l) = offset;
        end
        
        labels_discrete(i,:) = pos;
    %         d = (map_X-labels(i,1)).^2+(map_Y-labels(i,2)).^2;
    %         cur_gaussian = exp(-d/output_sigma/output_sigma)*1/sqrt(2*pi)/output_sigma;
    %         labels_gaussian(i,:)=cur_gaussian(:);
        if(mod(i,1000)==0)
            disp(i);
        end 
    end   
    stri = ['datasets/channel20_9x6_',num2str(j),'.mat'];
    save(stri,'channels','labels_discrete','offsets','cur_model','ap','-v7.3');
    clear
end
