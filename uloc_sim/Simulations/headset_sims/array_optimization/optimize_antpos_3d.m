%% Define global constants
addpath('../')
clearvars
%close all
pause('on')
theta_step = deg2rad(5);%0.01;
t_vals = (-deg2rad(180):theta_step:deg2rad(180)).'; % aoa values for FFT
p_vals = (-deg2rad(90):theta_step:deg2rad(90)).';
tp(:,1) = repmat(t_vals, length(p_vals), 1);
tmp = repmat(p_vals, 1,length(t_vals)).';
tp(:,2) = tmp(:);

% define opt
c = 3e8;
opt.freq = 4.4928e9;%6489.6e6;
opt.lambda = c./opt.freq;
opt.ant_sep = opt.lambda/2;
%radius = 2/8 * opt.lambda;

W = [cos(tp(:,2)).*cos(tp(:,1)), sin(tp(:,2)), cos(tp(:,2)).*sin(tp(:,1))];

% Provide some initial antenna position (or full array) ***************************
x = [0,-.0751202,0,0;0,0,-.0751202,0;0,0,0,-.0751202];

% Define search space SS
sdx = 0.005;
[ssx, ssy, ssz] = meshgrid(-.0751202:sdx:0, -.0751202:sdx:0, -.0751202:sdx:0); % Set size of search space
SS = [ssx(:),ssy(:),ssz(:)];
r_lowbound = 2/8.*opt.lambda; % Define minimum radius lowerbound for antenna separation
r_highbound = 4/8.*opt.lambda;

for iN = 1:6
    SSi = SS;
%     for in = 1:size(x,2)
%         SSi = [SSi; SS(vecnorm(SS-x(:,in).',2,2) <= r_highbound,:)];
%     end
    for in = 1:size(x,2)
        SSi = SSi(vecnorm(SSi-x(:,in).',2,2) >= r_lowbound,:);
    end
    ssi_len(iN) = length(SSi);
    Ax = W*SSi';
    vecs = exp(1j*2.*pi./opt.lambda.*Ax);

    X = exp(1j*2.*pi./opt.lambda.*W*x);
    Z = abs(sum(vecs'*X,2));
    [err_orth(iN), idx] = min(Z);
    
    x2 = SSi(idx,:).';
    x = [x, x2];
    
    % Evaluate Performance
%     X = exp(1j*2.*pi./opt.lambda.*W*x);
%     err(iN) = 0;
%     for ii = 1:length(tp)
%         den = 0;
%         for jj = setdiff(1:length(tp),ii)
%             den = den + abs(X(ii,:)*X(jj,:)')./size(x,2);
%         end
%         err(iN) = err(iN) + 1/den;
%     end
%     err(iN) = 1./(err(iN)./length(tp));
disp(iN)
end

%% Apply Kovarik's
if false
    X = exp(1j*2.*pi./opt.lambda.*W*x);
    Ak = X;
    I = eye(size(Ak,1));
    qk = 11; % Choose odd number...
    for ii = 1:10
        A = Ak*Ak';
        Kk = (I-A)*pinv(I+A);
        %K = zeros(size(Ak,1));
        %for jj = 0:qk
        %    K = K +(-1.*A).^jj;
        %end
        %Kk = (I-A)*K;
        Ak = (I+Kk)*Ak;
    end
end

%% Plots
figure
S = gen_auto_corr_steering(t_vals.', p_vals.', opt, x);

imagesc(S);

figure
scatter3(x(1,:),x(3,:),x(2,:))
xlabel('x')
ylabel('z')
zlabel('y')
axis equal

%save('arrays/test_ant_pos.mat', 'x', '-v7.3')

