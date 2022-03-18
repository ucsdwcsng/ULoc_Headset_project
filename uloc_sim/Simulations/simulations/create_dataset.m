clear;
% storage_path = './';
for j=2
    storage_path = '/media/user1/easystore/datasets/aditya_sim/';

    disp(j)
    file_name = [storage_path, sprintf('channels/channels_pix2pix_noerr_%d.mat',j)];
    %Load Channel file and initialize variables
%     load(file_name,'channels_wo_offset_noisy','channels_wo_offset', ...
%                 'channels_w_offset','labels','labels_noise', 'cur_model'); 
    load(file_name)
    
    input_grid_size=0.1;  %Grid size for input Gaussian
    output_grid_size=0.1;   %Grid size for output Gaussian
    output_sigma = 0.25; %Sigma for Gaussian output image
    theta_vals = -pi/2:0.01:pi/2; %Theta values to consider for multipath profiles
    d_vals = 0:0.1:30; %D vals to consider for distance profiles
    n_points=size(channels_w_offset, 1);
    n_lambda=length(cur_model.lambda);
    n_ap=length(ap);
    n_ant=length(ap{1});

    n_features=n_ap; % One distance profile per antenna, one angle profile per ap
    max_x =max(cur_model.walls(:,1));
    max_y =max(cur_model.walls(:,2));
    d1 = -max_x/2:input_grid_size:3/2*max_x;
    d2 = -max_y/2:input_grid_size:3/2*max_y;
    opt.freq = 3e8./(cur_model.lambda);
    opt.ant_sep = abs(ap{1}(2,1)-ap{1}(1,1));
    
    features_wo_offset = zeros(n_points,n_ap,length(d2),length(d1));
    features_w_offset = zeros(n_points,n_ap,length(d2),length(d1));
    
    features_wo_offset_noisy = zeros(n_points,n_ap,length(d2),length(d1));
    parfor i=1:n_points
        % Get features
        features_wo_offset(i,:,:,:) = generate_features_abs( ...
                        squeeze(channels_wo_offset(i,:,:,:)),ap, ...        
                                theta_vals,d_vals,d1,d2,opt);
        features_w_offset(i,:,:,:) = generate_features_abs( ...
                        squeeze(channels_w_offset(i,:,:,:)),ap, ...
                                theta_vals,d_vals,d1,d2,opt);
       
        if(mod(i,500)==0)
            disp(i);
        end    
    end
    labels_gaussian_2d=get_gaussian_labels_negative(labels,output_grid_size,0.25,d1,d2);

    disp(['Storing clean', j])
    stri = [storage_path, 'features/clean/dataset_pix_pix_',num2str(j),'.mat'];
    save(stri,'features_w_offset','features_wo_offset', ...
                'labels', ...
                'labels_gaussian_2d', ...
                'cur_model','ap','-v7.3');
    parfor i=1:n_points
        % Get noisy features
        features_wo_offset_noisy(i,:,:,:) = generate_features_abs( ...
                        squeeze(channels_wo_offset_noisy(i,:,:,:)),ap, ...
                                theta_vals,d_vals,d1,d2,opt);
       
        if(mod(i,500)==0)
            disp(i);
        end   
    end
    labels_gaussian_2d = get_gaussian_labels_negative(labels_noise,output_grid_size,0.25,d1,d2);

    features_wo_offset = features_wo_offset_noisy;
    labels = labels_noise;

    disp(['Storing noisy', j])
    stri = [storage_path, 'features/noisy/dataset_pix_pix_',num2str(j),'.mat'];
    save(stri,'features_w_offset','features_wo_offset', ...
                'labels', ...
                'labels_gaussian_2d', ...
                'cur_model','ap','-v7.3');

    clear
end
disp('done')