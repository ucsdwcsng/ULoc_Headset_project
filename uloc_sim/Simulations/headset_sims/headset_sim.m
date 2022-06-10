clearvars
%close all

root_folder = 'C:/users/tyler/Desktop/ULoc_Headset_project/';
addpath([root_folder,'uloc_sim/Simulations/simulations'])
addpath([root_folder,'uloc_sim/Simulations/legacy_sims'])

%% Define global constants

%pause('on')
theta_step = deg2rad(3);%0.01;
THETA_VALS = -deg2rad(90):theta_step:deg2rad(90); % aoa values for FFT
PHI_VALS = -deg2rad(90):theta_step:deg2rad(90);

% define opt
c = 3e8;
opt.cent_freq = 4.4928e9;%6489.6e6;
opt.freq = opt.cent_freq + 250e6*[-1:1/128:1-1/128];
opt.cent_lambda = c./opt.cent_freq;
opt.lambda = c./opt.freq;
%opt.ant_sep = 33.364e-3;
opt.ant_sep = opt.cent_lambda/2;
runtime = 0;

zeroFreq = mean(opt.freq);
% Time Resolution: If multipath arrives later than this time res relative 
% to direct path, multipath is ignored. 1e-9 for UWB, otherwise 1e9 to
% cover essentially any multipath
t_res = 1e9;
% TODO: weighted average gain patterns
%gain_pattern_az = load('/home/tyler/u_loc/az_gain_pattern.mat', 'gain_pattern');
%gain_pattern_e1 = load('/home/tyler/u_loc/az_gain_pattern.mat', 'gain_pattern');
gain_pattern_az.gain_pattern = ones(1,360);

%% Define the base_models
reflectors={};
walls = get_rectangle([-10,-10],20,20);  
ht = 10;
L_wall = 0;
R_wall = 5;

base_model.walls = walls; % just for aesthetics and plotting

base_model.reflectors = {};
base_model.obstacles = {};

base_model.lambda = opt.lambda;
base_model.amps = ones(length(base_model.lambda),1)*1000;
base_model.Ts = c/(80e6);
base_model.obs_attenuation = 1.0; % how much should the ray get attenuated when passing through an obstacle
num_sub_carriers = length(base_model.lambda);

%% Define Reflectors
% For reflectors, format is: reflectors{reflector set, reflector #} = [point1; point2; point3]
% Where points define a parallelogram plane in space, with point2 as the
% elbow of the 3 defined points.
num_sets_to_use = 0; % number of reflector sets to use (multipath)

%reflectors{1, 1} = [-2,3,1; 4,3,7; 4,7,7];
reflectors{1, 1} = [-2,0,2; -2,0,-2; 2,0,-2]; % Ground
%reflectors{1, 2} = [0.547, 1.933, 0.339; 0.547, 0.423, 0.339; 1.311,0.423,-0.545]; %RF2
wall = -0.2;
reflectors{1, 2} = [-0.5,0,wall ; 0.5,0,wall ; 0.5,3,wall];

%% Define the src and dst AP

% AP Locations
ap{1} = [0,0.5,0]; % x=-3 to offset

% Tag/Transmitter/Source Location
% body_space.mat has ground truth data storing possible tag locations
% around the human body.
load([root_folder,'uloc_sim/Simulations/headset_sims/body_space.mat']);
%src_center = ap{1}+[0.5,0,0]; 
%src = src_center;

n_init_ants = 4; % In the case of 3D options (5,6), n_init_ants is antennas per "axis". ie. 4 gives 10 total
arrtype = 3; % Select 5 for 3D-Axial, 6 for 3D-ZigZag, 9 to load in custom
% ArrTypes:
% Circular: 1, Coprime: 2, L: 3, Decahex: 4, 3D-Axial: 5 (n_init_ants refers to # antennas along each axis),
% 3D-ZigZag: 6, Optimized: 9 (From optimize_antpos_3d.m)
plt_env = true; % Plot Environment after simulation
plt_profile = true; % Plot AoA and Range per packet
ant_power = false; % Antenna power characterization (TODO)
recordVideo = false;
gen_path = true; % Use path from body_space.mat
localize = false; % Localize using range data
moving_head = false; % Headset rotates per packet or not (Need to add realistic rotations)
noise = true; % Add noise (Gaussian only currently, Multi-path disabled)
if recordVideo
    v = VideoWriter('Anim_Frames.avi'); v.FrameRate = 1; open(v);
end
% ************ Mostly don't need to edit after this point...************************************************************************************************

% For now, only using azimuth (uniform) gain
aoa_err = [];

src_orient = 90;

dst_orient_az = -90;
if dst_orient_az < 0
    dst_orient_az = dst_orient_az + 360;
end

TX_ANT_NUM = 1;%size(src, 1);
%{
dst_orient_e1 = rad2deg(atan2(dst_e1(4,2)-dst_e1(1,2), dst_e1(4,3)-dst_e1(1,3)))-90;
if dst_orient_e1 < 0
    dst_orient_e1 = dst_orient_e1 + 360;
end
%}

% Generate a path for the src, or default to wherever src was placed
if gen_path
     path = gt_data;
     path = path(1:500,:);
     if moving_head
         ap_locs = [0.5*(rand(length(path),1)-0.5)+path(:,1),...
                    mean(path(:,2))*ones(length(path),1),...
                    0.5*rand(length(path),1)+path(:,3)];
     end
else
    path = src;
    if moving_head
        ap_locs(1,:) = ap{1};
    end
end
totallen = size(path,1);

if ~moving_head
    rots = 0;
    array = gen_array(n_init_ants, arrtype, opt, ap, rots);
    n_ants = size(array.local_ants{1,1},1);
    R = vecnorm(path-ap{1},2,2);
end

if noise % Initialize Random Multipath Noise (TODO: Add across CIR taps)
    array = gen_array(n_init_ants, arrtype, opt, ap, 0);
    const = 1j*2*pi./opt.cent_lambda;
    ant_pos = array.local_ants{1}.';
    N = [];
    for jj=1:length(PHI_VALS)
        th = THETA_VALS;
        phi = PHI_VALS(jj);
        wave_vec = const * [cos(phi)*cos(th); sin(phi)*ones(1,length(THETA_VALS)); cos(phi)*sin(th)];            
        steering_vec = exp(wave_vec.' * -ant_pos);%exp(dot(wave_vec, ant_pos)).';
        N(1+(jj-1)*size(THETA_VALS,2):jj*size(THETA_VALS,2),:) = steering_vec;
    end
end

for it = 1:totallen % Main loop
%% Define the channels
if moving_head % Add head rotations
    ap{1} = ap_locs(it,:);
    rots = randn;
    array = gen_array(n_init_ants, arrtype, opt, ap, rots);
    n_ants = size(array.local_ants{1,1},1);
    R = vecnorm(path-ap{1},2,2);
end

% Simulate Channel
channels{1} = zeros(num_sub_carriers, array.RX_ANT_NUM, TX_ANT_NUM);
rays{1} = {};
are_rays_blocked{1} = {};
for rx_ant=1:array.RX_ANT_NUM
    for tx_ant = 1:TX_ANT_NUM
        model = base_model;
        for ii=1:size(reflectors, 2)
            if (~isempty(reflectors))
                model.reflectors = [model.reflectors,  {reflectors{1:num_sets_to_use, ii}}];
            else
                model.reflectors =  {reflectors{1:num_sets_to_use, ii}};
            end
        end
        [channels{1}(:, rx_ant, tx_ant), r, blocked] = ...
                get_channels_from_model(model, path(it, :), array.ants{1}(rx_ant, :), t_res, false, src_orient, dst_orient_az, gain_pattern_az.gain_pattern); 
        rays{1,1}{rx_ant, tx_ant} = r;
        are_rays_blocked{1,1}{rx_ant, tx_ant} = blocked; 
    end
end

%% Theta-Phi Pred
for iap = 1:length(ap) % Only 1 AP currently, loop unecessary
    H = channels{iap}.';
    H = H./abs(H); % Normalize antenna power
    if noise
        H = H + 0.75*(randn(n_ants,1) + 1j*randn(n_ants,1)); % Gaussian Noise
        H = H./abs(H);
        ii = 0;% Number of random paths. If 0, no random MP added
        for jj = 1:ii
            noise_ants = rand(n_ants,1);
            tmp = sort(noise_ants);
            noise_ants = noise_ants.*((noise_ants - tmp(5)) <= 0);
            noise_ch = noise_ants .* N(ceil(rand(1)*size(N,1)),:).';
            %H = H + ( 0.5 + 0.5*rand(1) )*noise_ch;
            H = H + 0.75*noise_ch;
        end
        H = H./max(abs(H));
    end
    
    % Generate AoA Profile
    %[P,A,dt] = gen_theta_phi_fft_general(H, THETA_VALS, PHI_VALS, opt, array.local_ants{iap}.',plt_profile);
    [P,dt] = gen_theta_phi_music_general(H, THETA_VALS, PHI_VALS, opt, array.local_ants{iap}.',1,plt_profile);
    runtime = runtime + dt;
    [tmp, tmax] = max(abs(P));
    [~, pmax] = max(tmp);
    aoa_pred = [rad2deg(THETA_VALS(tmax(pmax))), rad2deg(PHI_VALS(pmax))];
    
    % Calculate Ground Truth AoA
    if arrtype < 5 % 2D Cases
        gt_aoa(iap,:) = rad2deg( [atan2(path(it,3)-ap{iap}(3), path(it,1)-ap{iap}(1)),...
                                 atan2(path(it,2)-ap{iap}(2), norm(path(it,[1,3]) - ap{iap}([1,3])))] );
        gt_aoa(iap,1) = gt_aoa(iap,1) - rad2deg(rots(iap));

        if gt_aoa(iap, 1) > 90
            gt_aoa(iap, 1) = gt_aoa(iap,1) - 180;
        elseif gt_aoa(iap, 1) < -90
            gt_aoa(iap, 1) = gt_aoa(iap,1) + 180;
        end
    else % 3D Cases
        gt_aoa(iap,:) = rad2deg( [atan2(path(it,3)-ap{iap}(3), path(it,1)-ap{iap}(1)),...
                                 atan2(path(it,2)-ap{iap}(2), norm(path(it,[1,3]) - ap{iap}([1,3])))] );
        gt_aoa(iap,1) = gt_aoa(iap,1) + rad2deg(rots(iap));
        
        if gt_aoa(iap, 1) > 180
            gt_aoa(iap, 1) = gt_aoa(iap,1) - 360;
        elseif gt_aoa(iap, 1) < -180
            gt_aoa(iap, 1) = gt_aoa(iap,1) + 360;
        end
    end
    
    if plt_profile % Plot Profile
        figure(1)
        subplot(1,2,1)
        theta_phi_plot = pcolor(rad2deg(PHI_VALS),rad2deg(THETA_VALS),abs(P));
        set(theta_phi_plot, 'EdgeColor', 'none');
        hold on
        xline(aoa_pred(2),'r')
        yline(aoa_pred(1),'r')
        xline(gt_aoa(iap,2),'k')
        yline(gt_aoa(iap,1),'k')
        hold off
        xlabel('Phi')
        ylabel('Theta')
        title(['AP: ',num2str(iap), ' Err: ', num2str(gt_aoa(iap,1)-aoa_pred(1)), ', ', num2str(gt_aoa(iap,2)-aoa_pred(2)), ' n antennas: ',num2str(n_ants)])

%         if ant_power
%             figure(2)
%             theta_phi_plot = pcolor(rad2deg(PHI_VALS),rad2deg(THETA_VALS),abs(A));
%             set(theta_phi_plot, 'EdgeColor', 'none');
%             hold on
%             xline(aoa_pred(2),'r')
%             yline(aoa_pred(1),'r')
%             xline(gt_aoa(iap,2),'k')
%             yline(gt_aoa(iap,1),'k')
%             hold off
%             xlabel('Phi')
%             ylabel('Theta')
%             colorbar
%         end
        waitforbuttonpress
    end
    aoa_err = [aoa_err; gt_aoa(iap,:) - aoa_pred];
end

%% Localization
if localize % Need to change
    aoa_extrapolation(iap,:) = [R(it)*cos(deg2rad(aoa_pred(2)))*cos(deg2rad(aoa_pred(1))),...
                            R(it)*sin(deg2rad(aoa_pred(2))),...
                            R(it)*cos(deg2rad(aoa_pred(2)))*sin(deg2rad(aoa_pred(1)))];

    aoa_extrapolation(iap,:) = aoa_extrapolation(iap,:)*inv(squeeze(array.M(iap,:,:))) + ap{iap};
    loc_pred(it,:) = aoa_extrapolation(iap,:);
    loc_err(it) = norm(loc_pred(it,:)-path(it,:));
end
if mod(it,floor(totallen/10)) == 0
    disp([num2str(it),'/',num2str(totallen)]);
end

%sgtitle(['Loc Err: ',num2str(loc_err(it))])
if recordVideo
    frame = getframe(gcf);
    writeVideo(v, frame);
end
%pause(0.75)

% Record Metrics
metrics.AoA_err(it,:) = aoa_err(it,:);
metrics.AoA_err_norm(it) = norm(aoa_err(it,:));
[~, img, peaks] = find_peaks(P, THETA_VALS, PHI_VALS);
bw = imbinarize(img,'global');
props = regionprops(bw, 'Centroid','Area','Eccentricity');
p=[]; pk=[];
for ik = 1:length(props)
    p(ik) = norm(rad2deg([THETA_VALS(round(props(ik).Centroid(2))),PHI_VALS(round(props(ik).Centroid(1)))]) - gt_aoa);
end
for ik = 1:length(peaks{1})
    pk(ik) = norm([peaks{1}(ik), peaks{2}(ik)] - gt_aoa);
end
[~,idx] = min(p);
[~,idx2] = min(pk);

metrics.close_peak(it,:) = [peaks{1}(idx2), peaks{2}(idx2)];%rad2deg([THETA_VALS(round(props(idx).Centroid(2))),PHI_VALS(round(props(idx).Centroid(1)))]);
metrics.num_peaks(it) = length(props);
metrics.peak_size(it) = props(idx).Area;
metrics.eccentricity(it) = props(idx).Eccentricity;
metrics.close_peak_err(it,:) = (gt_aoa - metrics.close_peak(it,:));
metrics.close_peak_err_norm(it) = norm(metrics.close_peak_err(it,:));

if length(props)>1
    for ik = 1:length(props)
        p(ik) = norm([peaks{1}(ik), peaks{2}(ik)] - metrics.close_peak(it,:));
    end
    p=p(p~=0);
    metrics.min_dist(it) = min(p);
else
    metrics.min_dist(it) = 0;
end
end % End Loop for 1 iteration
runtime = runtime ./totallen;
if recordVideo
    close(v);
end

%% Plot Environment
if plt_env
    figure
    hold on
    scatter(path(:,3),path(:,1), 100,'k', '.')
    if localize
        scatter(loc_pred(:,3),loc_pred(:,1), 100, loc_err,'.')
    end
    for iap = 1:length(ap)
        scatter(array.ants{iap}(:,3),array.ants{iap}(:,1))
    end
    if localize
        legend('Path', 'Pred','AP1')
    else
        legend('Path','AP1')
    end
    if num_sets_to_use
    for iref = 1:size(reflectors,2)
        plot(reflectors{1,iref}(:,3), reflectors{1,iref}(:,1),'-*','color','k')
    end
    end
    hold off
    %title('Optimized AP Error Along Path')
    title('Optimized AP Setup')
    colorbar
    axis equal
    xlabel('Z'); ylabel('X');
end
aoa_err = abs(aoa_err);
disp(['Median AoA Err(t,p): ', num2str(prctile(aoa_err(:,1), 50)), ', ', num2str(prctile(aoa_err(:,2), 50))])
disp(['90th AoA Err(t,p): ', num2str(prctile(aoa_err(:,1), 90)), ', ', num2str(prctile(aoa_err(:,2), 90))])
disp(['99th AoA Err(t,p): ', num2str(prctile(aoa_err(:,1), 99)), ', ', num2str(prctile(aoa_err(:,2), 99))])

if localize
    disp(['Median XYZ Err: ', num2str(prctile(loc_err, 50))])
    disp(['90th XYZ Err: ', num2str(prctile(loc_err, 90))])
    disp(['99th XYZ Err: ', num2str(prctile(loc_err, 99))])
end

%S = gen_auto_corr_steering(THETA_VALS, PHI_VALS, opt, array.local_ants{1}.');
