clear
close all
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
walls = get_rectangle([0,0],20,15);  % Walls are defined just for plotting. They are non-reflecting, non-blocking                              
obstacles{1} = get_rectangle([7,5],2,3); %Obstacles are bounded and opaque to wireless signals
%obstacles{2} = get_rectangle([16,5],2,3); %get rectangle([x,y],w,h) returns a rectangle at x,y with width w and height h
reflectors{1} = [20,2;20,13]; % Reflectors are linear for ease 
model.walls = walls;
model.reflectors =  reflectors;
model.obstacles = obstacles;
model.obstacles={};
model.reflectors={};
model.lambda = 3e8./(2.4e9:1e6:2.48e9);
model.amps = ones(length(model.lambda),1)*100; % Lambda and amps can be arrays
%% Define the setup
max_x = max(model.walls(:,1));
max_y = max(model.walls(:,2));
ant_sep = min(model.lambda)/2;
opt.freq=2.4e9:1e6:2.48e9;
opt.ant_sep=ant_sep;
ap{1}=[10,0;
    10+ant_sep,0;
    10+2*ant_sep,0;
    10+3*ant_sep,0;
    10+4*ant_sep,0];

ap{2}=[10,15;
    10+ant_sep,15;
    10+2*ant_sep,15;
    10+3*ant_sep,15;
    10+4*ant_sep,15];

ap{3}=[0,7.5-2*ant_sep;
        0,7.5-ant_sep;
        0, 7.5;
        0,7.5+ant_sep;
        0,7.5+2*ant_sep];

ap{4}=[20,7.5-2*ant_sep;
        20,7.5-ant_sep;
        20, 7.5;
        20,7.5+ant_sep;
        20,7.5+2*ant_sep];
ap_orient=[1,1,2,2];
n_total = 1;
n_ant_per_ap = 5;
theta_vals = [-pi/2:0.01:pi/2];
d_vals = -15:1:15;
d1 = 0:0.5:max_x;
d2=0:0.5:max_y;
%% Generate the points. Ap 1 antenna 1 serves as the reference
%features=zeros(n_total,length(model.lambda)*length(ap)*n_ant_per_ap*2);
features = zeros(n_total,length(ap),length(ap),length(d2),length(d1));
labels = zeros(n_total,2);
channels = zeros(length(model.lambda),length(ap),n_ant_per_ap);
channels_rel= zeros(length(model.lambda),length(ap),n_ant_per_ap,length(ap)-1);
rng;
plotthis = 1;
start_time = now;
%parpool(4);
for i=1:n_total
    pos = [rand()*max(model.walls(:,1)),rand()*max(model.walls(:,2))];
%     pos = [10, 7.5];
    for j=1:length(ap)
        for k=1:n_ant_per_ap
            channels(:,j,k)=get_channels_from_model(model,pos,ap{j}(k,:),false);
            ind = 1:1:length(ap);
            ind(j) = [];
            for l = 1:length(ap)-1
                channels_rel(:,j,k,l) = get_channels_from_model(model,pos,ap{j}(k,:),false).*conj(get_channels_from_model(model,pos,ap{ind(l)}(1,:),false));
                channels_rel(:,j,k,l) = awgn(squeeze(channels_rel(:,j,k,l)),30);
            end
        end
        
        channels(:,j,:) = awgn(squeeze(channels(:,j,:)),30);
        for l = 1:length(ap)-1
            
            Prel =  compute_multipath_profile2d_fast_edit(squeeze(channels_rel(:,j,:,l)),theta_vals,d_vals,opt);
            Prel_out = convert_relative_spotfi_to_2d_edit(Prel,ap{j},ap{ind(l)},theta_vals,d_vals,d1,d2);
            Prel2 =  compute_multipath_profile2d_fast(squeeze(channels_rel(:,j,:,l)),theta_vals,d_vals,opt);
            Prel_out2 = convert_relative_spotfi_to_2d_edit(Prel2,ap{j},ap{ind(l)},theta_vals,d_vals,d1,d2);
            figure(j),subplot(2,2,l),imagesc(db(abs(Prel_out))),title(['Ap ' num2str(j) ' with respective to Ap ' num2str(ind(l))]);
            figure(j+length(ap)),subplot(2,2,l),imagesc(db(abs(Prel_out2))),title(['Ap ' num2str(j) ' with respective to Ap ' num2str(ind(l))]);
            features(i,j,l,:,:) = abs(Prel_out);
            
        end
        P = compute_multipath_profile2d_fast(squeeze(channels(:,j,:)),theta_vals,d_vals,opt);
        P_out = convert_spotfi_to_2d(P,theta_vals,d_vals,d1,d2,ap{j});
        P1 = compute_multipath_profile2d_fast_edit(squeeze(channels(:,j,:)),theta_vals,d_vals,opt);
        P_out1 = convert_spotfi_to_2d(P1,theta_vals,d_vals,d1,d2,ap{j});
        %channels = channels./repmat(channels(:,1,1),1,length(ap),n_ant_per_ap);
        figure(j),subplot(2,2,4),imagesc(db(abs(P_out1))),title(['Ap ' num2str(j)]);
        figure(j+length(ap)),subplot(2,2,4),imagesc(db(abs(P_out))),title(['Ap ' num2str(j)]);
        features(i,j,length(ap),:,:) = abs(P_out);
        plotthis = plotthis.*P_out;
    end
    figure,imagesc(d1,d2,db(plotthis));
    labels(i,:) = pos;
    if(mod(i,100)==0)
        disp([i,(now -start_time)*24*60]);
    end
end

%delete(gcp('nocreate'));