clearvars
close all
root_folder = 'C:/users/tyler/Desktop/ULoc_Headset_project/';
addpath([root_folder,'uloc_sim/Simulations/simulations'])
addpath([root_folder,'uloc_sim/Simulations/legacy_sims'])

use_matlab_imu = false;
gen_new_traj = false;
g = 9.80665;

%% Generate Trajectory
if gen_new_traj
    imu_Fs = 2000;
    len = 2000;
    init_xyz = [cos((1:len)*.01); sin((1:len)*0.1); (-len/2:len/2-1)*.001];
    trajectory = waypointTrajectory(init_xyz.', 0.02*(1:len));
    trajectory.SampleRate = imu_Fs;%1./refresh_rate;

    %% Initialize IMU
    if use_matlab_imu
        imuFs = 1000; % Up to 8000Hz
        imu = imuSensor('accel-gyro-mag', 'SampleRate', imuFs);
        %imu.MagneticField = [19.5281 -5.0741 48.0067];

        % Accelerometer
        imu.Accelerometer.MeasurementRange =  100;%2*g;
        imu.Accelerometer.Resolution = (g./16384); %0.0023928;
        imu.Accelerometer.ConstantBias = 0;%0.19;
        imu.Accelerometer.NoiseDensity = 400./(g*1000);%0.0012356;

        % Gyroscope
        imu.Gyroscope.MeasurementRange = deg2rad(250);
        imu.Gyroscope.Resolution = deg2rad(1./13100); %deg2rad(0.0625);
        imu.Gyroscope.ConstantBias = deg2rad(0);
        imu.Gyroscope.AxesMisalignment = 0;%2;
        imu.Gyroscope.NoiseDensity = deg2rad(0.000005);%.005 %deg2rad(0.025);

        % Magnetometer: MPU6050 uses an unknown 3rd party Magnetometer, use defaults
        imu.Magnetometer.MeasurementRange = 1000;
        imu.Magnetometer.Resolution = 0.1;
        imu.Magnetometer.ConstantBias = 0;%100;
        imu.Magnetometer.NoiseDensity = 0.3/ sqrt(50);
    end
    %% Generate IMU Readings
    ii = 1;
    while ~isDone(trajectory)% = 1:len
        [pos,orient,vel,acc,angvel] = trajectory();
        if ii == 1
            tpos = []; torient = quaternion(); tvel = []; tacc = []; tangvel = [];
        else
            tpos = [tpos;pos];
            torient = [torient;orient];
            tvel = [tvel;vel];
            tacc = [tacc;acc];
            tangvel = [tangvel;angvel];
        end
        ii = ii+1;
    end
    % tacc = tacc./1000; % Returns in mm/s^2 for some reason
    % tacc = tacc./refresh_rate; % tacc is currently just diff of velocity, needs additional division
    % tangvel = deg2rad(tangvel); % Returns as deg for some reason
    if use_matlab_imu
        [imu_acc, imu_gyro, imu_mag] = imu(tacc,tangvel,torient);
        gt_pos = tpos; gt_orient = quat2rotm(torient); gt_vel = tvel; gt_acc = tacc; gt_angvel = tangvel;
    else
        torient = quat2rotm(torient);
        for ii = 1:size(tacc,1)
            imu_acc(ii,:) = tacc(ii,:)*squeeze(torient(:,:,ii));
            imu_gyro(ii,:) = tangvel(ii,:)*squeeze(torient(:,:,ii));
        end
        gt_pos = tpos; gt_orient = torient; gt_vel = tvel; gt_acc = tacc; gt_angvel = tangvel;
    end

    t = 1./imu_Fs * (1:length(gt_pos));
    refresh_rate = 1./imu_Fs;
    
    % AoA and Range Ground Truth
    origin = [0,0,0];
    for ii = 1:length(t)
        gt_aoa(ii,:) = [ atan2(gt_pos(ii,2),gt_pos(ii,1)) ...
                   atan2(gt_pos(ii,3), norm(gt_pos(ii,[1,2])) )]; % Theta, Phi
        gt_range(ii) = norm(gt_pos(ii,:)-origin);
    end
    save('Sim_IMU_Traj.mat','gt_pos','gt_orient','gt_vel','gt_acc','gt_angvel','imu_acc','imu_gyro','imu_Fs','t','refresh_rate','gt_aoa','gt_range','-v7.3');
else
    load('Sim_IMU_Traj.mat')
end
imu_gyro = imu_gyro + deg2rad(0.05) * randn(size(imu_gyro)); % See MPU6050 Datasheet
imu_acc = imu_acc + 0.0126491106407 * randn(size(imu_acc)); % See MPU6050 Datasheet, 400./1e6 *sqrt(1000)

range_est = gt_range + 0.02*randn(size(gt_range));
max_aoa_err = deg2rad( 5 );
aoa_est = gt_aoa + 2*max_aoa_err*rand(size(gt_aoa))-max_aoa_err;

%% AoA Analysis Without IMU
for ierr = 0:40
    range_est = gt_range + 0.02*randn(size(gt_range));
    max_aoa_err = deg2rad( ierr );
    aoa_est = gt_aoa + 2*max_aoa_err*rand(size(gt_aoa))-max_aoa_err;
    for ii = 1:size(aoa_est,1)
        xyz_pred(ii,:) = [cos(aoa_est(ii,2))*cos(aoa_est(ii,1)),...
                    cos(aoa_est(ii,2))*sin(aoa_est(ii,1)),...
                    sin(aoa_est(ii,2))] * range_est(ii);
    end
    err(ierr+1,:) = vecnorm(gt_pos-xyz_pred,2,2); % Loc Err
end

% Plot
figure
for ii = [1,6,11,16,21,26,31,36,41]
    hold on
    ecdf(err(ii,:))
end
grid on
grid minor
xlim([0,0.5])
xlabel('Error(m)')
xlabel('Error (m)')
title('CDF Localization with Noisy Range(2cm std) and Increasing AoA Error')
legend('0\circ','5\circ','10\circ','15\circ','20\circ','25\circ','30\circ','35\circ','40\circ')

%% Kalman Filter (IMU Only)
% Define:
%   y: Unfiltered estimate of pose
%   x: State space Matrix involving: [x0;v0;acc0] 3x3

Q=1; % Process Noise Covariance
R=0.05; % Sensor Noise Covariance
A = eye(3); % State to State
B = [1;1;1]; % Input to State
C = [1 refresh_rate 0.5*refresh_rate.^2]; % State to Output
D=0; % Input to Output (Feedthrough)

P = B*Q*B'; % Initial Error Covariance
x = zeros(3); % Initial State Condition

y_est = zeros(length(t),3); % Filtered Response
y = zeros(length(t),3); % Naive Response
errcov = zeros(length(t),1); % Error Covariance, should reach a steady state over time

% Initialize filter
%FUSE = imufilter(); FUSE.OrientationFormat = 'Rotation matrix';
%FUSE.AccelerometerNoise=1e-13;FUSE.GyroscopeNoise=1e-13;FUSE.GyroscopeDriftNoise=1e-13;FUSE.LinearAccelerationNoise=1e-13;Fuse.LinearAccelerationDecayFactor=0;
start = 39;
x = [gt_pos(start,:);gt_vel(start,:);gt_acc(start,:)];
for ii = start:length(t)
    % Naive Estimate of y
    if ii == start
        orient = squeeze(gt_orient(:,:,ii)); % Assume know initial orientation
        orient_sv(:,:,ii) = orient;
        vi = gt_vel(ii,:); % Assume know initial velocity
        v_naive(ii,:) = vi;
        acc_sv(ii,:) = imu_acc(ii,:)*pinv(orient);
        y(ii,:) = gt_pos(ii,:);% + 0.5*(imu_acc(ii,:)*pinv(orient))*refresh_rate.^2 + vi*refresh_rate; % Assume know inital position
    else
        y(ii,:) = y_est(ii-1,:) + (0.5*(imu_acc(ii-1,:)*pinv(orient))*refresh_rate.^2) + vi*refresh_rate;
        
        vi = vi + (imu_acc(ii-1,:)*pinv(orient))*refresh_rate;
        v_naive(ii,:) = vi;
        angs = imu_gyro(ii-1,:)*refresh_rate;
        orient_z = [cos(angs(3)) -sin(angs(3)) 0; sin(angs(3)) cos(angs(3)) 0; 0 0 1];
        orient_y = [cos(angs(2)) 0 sin(angs(2)); 0 1 0; -sin(angs(2)) 0 cos(angs(2))];
        orient_x = [1 0 0; 0 cos(angs(1)) -sin(angs(1)); 0 sin(angs(1)) cos(angs(1))];     
        orient = orient_z * orient_y * orient_x * orient;

        orient_sv(:,:,ii) = orient;
        acc_sv(ii,:) = imu_acc(ii,:)*pinv(orient);
    end
    % Measurement Update
    Mxn = P*C'/(C*P*C'+R);        % Mx[n]
    x = x + Mxn*(y(ii,:)-C*x);      % x[n|n]
    P = (eye(length(x))-Mxn*C)*P; % P[n|n]
    
    y_est(ii,:) = C*x;
    errcov(ii) = C*P*C';
    
    % Time Update
    %x = A*x + B*u(ii);   % x[n+1|n]
    x = [y_est(ii,:); vi; (0.5*imu_acc(ii,:)*refresh_rate.^2)*pinv(orient)];
    P = A*P*A' + B*Q*B'; % P[n+1|n]
end

%% Plot
figure
subplot(3,1,1)
plot(t,gt_pos(:,1),'--')
hold on
plot(t,y_est(:,1))
grid on
grid minor
xlabel('Time(s)')
ylabel('X')
legend('Ground Truth', 'IMU Integrated')
subplot(3,1,2)
plot(t,gt_pos(:,2),'--')
hold on
plot(t,y_est(:,2))
grid on
grid minor
xlabel('Time(s)')
ylabel('Y')
subplot(3,1,3)
plot(t,gt_pos(:,3),'--')
hold on
plot(t,y_est(:,3))
grid on
grid minor
xlabel('Time(s)')
ylabel('Z')
sgtitle('Noisy XYZ from IMU only, KF')

%% Kalman Filter, AoA/Range and IMU
aoa_errs = 10;%[0:5:40];
for ierr = 1:length(aoa_errs)
    range_est = gt_range + 0.02*randn(size(gt_range));
    max_aoa_err = deg2rad( aoa_errs(ierr) );
    aoa_est = gt_aoa + 2*max_aoa_err*rand(size(gt_aoa))-max_aoa_err;
    
    % Define:
    %   y: Unfiltered estimate of pose
    %   x: State space Matrix involving: [x0;v0;acc0;x_uloc] 4x3
    Q=diag([2.765, 0.0028, 0.0481, cov(err(aoa_errs(ierr)+1,39:end))]);%1; % Process Noise Covariance
    R=0; % Sensor Noise Covariance
    A = eye(4); % State to State
    B = [1;1;1;1]; % Input to State
    C = [1 refresh_rate 0.5*refresh_rate.^2 1]; % State to Output
    D=0; % Input to Output (Feedthrough)

    P = Q;%B*Q*conj(B); % Initial Error Covariance
    x = zeros(4,3); % Initial State Condition

    y_est = zeros(length(t),3); % Filtered Response
    y = zeros(length(t),3); % Naive Response
    errcov = zeros(length(t),1); % Error Covariance, should reach a steady state over time

    % Initialize filter
    %FUSE = imufilter(); FUSE.OrientationFormat = 'Rotation matrix';
    %FUSE.AccelerometerNoise=1e-13;FUSE.GyroscopeNoise=1e-13;FUSE.GyroscopeDriftNoise=1e-13;FUSE.LinearAccelerationNoise=1e-13;Fuse.LinearAccelerationDecayFactor=0;
    start = 39;

    uloc_xyz(start,:) = [cos(aoa_est(start,2))*cos(aoa_est(start,1)),...
                        cos(aoa_est(start,2))*sin(aoa_est(start,1)),...
                        sin(aoa_est(start,2))] * range_est(start);
    x = [gt_pos(start,:);gt_vel(start,:);gt_acc(start,:);uloc_xyz(start,:)];
    for ii = start:length(t)
        % Naive Estimate of y
        if ii == start
            orient = squeeze(gt_orient(:,:,ii)); % Assume know initial orientation
            orient_sv(:,:,ii) = orient;
            vi = gt_vel(ii,:); % Assume know initial velocity
            v_naive(ii,:) = vi;
            acc_sv(ii,:) = imu_acc(ii,:)*pinv(orient);

            y(ii,:) = gt_pos(ii,:);% + 0.5*(imu_acc(ii,:)*pinv(orient))*refresh_rate.^2 + vi*refresh_rate; % Assume know inital position
        else
            y(ii,:) = y_est(ii-1,:) + (0.5*(imu_acc(ii-1,:)*pinv(orient))*refresh_rate.^2) + vi*refresh_rate;

            vi = vi + (imu_acc(ii-1,:)*pinv(orient))*refresh_rate;
            v_naive(ii,:) = vi;
            angs = imu_gyro(ii-1,:)*refresh_rate;
            orient_z = [cos(angs(3)) -sin(angs(3)) 0; sin(angs(3)) cos(angs(3)) 0; 0 0 1];
            orient_y = [cos(angs(2)) 0 sin(angs(2)); 0 1 0; -sin(angs(2)) 0 cos(angs(2))];
            orient_x = [1 0 0; 0 cos(angs(1)) -sin(angs(1)); 0 sin(angs(1)) cos(angs(1))];
            orient = orient_z * orient_y * orient_x * orient;

            orient_sv(:,:,ii) = orient;
            acc_sv(ii,:) = imu_acc(ii,:)*pinv(orient);
        end
        % Measurement Update
        Mxn = P*C'/(C*P*C'+R);        % Mx[n]
        x = x + Mxn*(y(ii,:)-(C*x)./2);      % x[n|n]
        P = (eye(length(x))-Mxn*C)*P; % P[n|n]

        y_est(ii,:) = (C*x)./2;
        errcov(ii) = C*P*C';

        % Time Update
        %x = A*x + B*u(ii);   % x[n+1|n]
        uloc_xyz(ii,:) = [cos(aoa_est(ii,2))*cos(aoa_est(ii,1)),...
                        cos(aoa_est(ii,2))*sin(aoa_est(ii,1)),...
                        sin(aoa_est(ii,2))] * range_est(ii);
        x = [y_est(ii,:);...
             vi;...
             (0.5*imu_acc(ii,:)*refresh_rate.^2)*pinv(orient);...
             uloc_xyz(ii,:)];
        P = A*P*A' + Q;%B*Q*conj(B); % P[n+1|n]
    end
    loc_err(ierr,:) = vecnorm(gt_pos-y_est,2,2);
end

% Without KF, combine IMU and ULoc (Simple Averaging)
for ii = 39:length(t)
    if ii==39
        av_est(ii,:) = gt_pos(ii,:);
    else
        av_est(ii,:) = (av_est(ii-1,:)+(refresh_rate*v_naive(ii-1,:)+0.5*acc_sv(ii-1,:)*refresh_rate.^2) +...
                       uloc_xyz(ii,:))./2;
    end
end
naive_av_err = vecnorm((av_est-gt_pos),2,2);

%% Plot
figure
subplot(3,2,1)
plot(t,gt_pos(:,1),'--')
hold on
plot(t,y_est(:,1))
grid on
grid minor
xlabel('Time(s)')
ylabel('X')
legend('Ground Truth', 'AoA, IMU')
subplot(3,2,3)
plot(t,gt_pos(:,2),'--')
hold on
plot(t,y_est(:,2))
grid on
grid minor
xlabel('Time(s)')
ylabel('Y')
subplot(3,2,5)
plot(t,gt_pos(:,3),'--')
hold on
plot(t,y_est(:,3))
grid on
grid minor
xlabel('Time(s)')
ylabel('Z')
subplot(1,2,2)
ecdf(loc_err(39:end))
hold on
ecdf(naive_av_err(39:end))
ecdf(err(11,39:end))
hold off
grid minor
grid on
legend('KF', 'No KF', 'ULoc')
title('XYZ Error')
sgtitle('Noisy IMU, AoA, KF')