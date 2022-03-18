%% Define global constants
clearvars
%close all
pause('on')
theta_step = deg2rad(0.5);%0.01;
t_vals = (-deg2rad(90):theta_step:deg2rad(90)).'; % aoa values for FFT
p_vals = 0;
tp(:,1) = repmat(t_vals, length(p_vals), 1);
tp(:,2) = repmat(p_vals, length(t_vals), 1);

% define opt
c = 3e8;
opt.freq = 4.4928e9;%6489.6e6;
opt.lambda = c./opt.freq;
opt.ant_sep = opt.lambda/2;

W = [sin(tp(:,1)),cos(tp(:,1))];

% Provide some initial antenna position (or full array) ***************************
%x0 = [0,0.5342;0,0];

%
sep = opt.lambda*4/8;
x0 = unique([3*sep*[0:4], 4*sep*[0:3]]);
x0 = [x0;zeros(1,length(x0))];
%}
x = [0,x0(1,end);0,0];
x = [0;0];
% Define search space SS
%box_len = [-0.095:0.005:0.1];
sdx = 0.005;
%[ssx, ssy] = meshgrid(0:sdx:x0(1,end), 0);
[ssx, ssy] = meshgrid(0:sdx:0.125, 0);
SS = [ssx(:),ssy(:)];
r_lowbound = 0/8.*opt.lambda; % Define minimum radius lowerbound for antenna separation
r_highbound = 4/8.*opt.lambda;

for iN = 1:5
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
    
    %x2 = real( Winv * (log(vecs(:,idx))./(1j*2.*pi./opt.lambda)) );
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

%% Plots
% figure
% S = gen_auto_corr_steering(t_vals.', p_vals.', opt, x);
% imagesc(S);
X = exp(1j*2.*pi./opt.lambda.*W*x);
X0 = exp(1j*2.*pi./opt.lambda.*W*x0);
AX = abs(X*X');
AX0 = abs(X0*X0');
on_diag1 = mean(diag(AX));
on_diag2 = mean(diag(AX0));
score(1) = sum(on_diag1./(sum(AX,2)-on_diag1))./length(AX);
score(2) = sum(on_diag2./(sum(AX0,2)-on_diag2))./length(AX0);

figure
subplot(2,2,1);imagesc(tp(:,1),tp(:,1),AX);title(['Optimized 8 Antenna, Score: ',num2str(score(1))])
subplot(2,2,3);imagesc(tp(:,1),tp(:,1),AX0);title(['Literature 8 Antenna, Score: ',num2str(score(2))])


subplot(1,2,2)
scatter(x(1,:),x(2,:))
hold on
scatter(x0(1,:),x0(2,:),'r*')
legend('Optimized','Literature')
xlabel('x')
ylabel('z')
zlabel('y')
axis equal
%save('test_ant_pos_1d_opt.mat', 'x', '-v7.3')

%% Calculate DoF
xd0 = x0(1,:);
xd1 = x(1,:);
d0 = [];
d1 = [];
for ii = 1:length(xd0)
    d0 = [d0,xd0-xd0(ii)];
    d1 = [d1,xd1-xd1(ii)];
end
d0 = unique(round(8*unique(abs(d0))./(opt.lambda*1/2))./8);
d1 = unique(round(8*unique(abs(d1))./(opt.lambda*1/2))./8);
d0 = d0(2:end);
d1 = d1(2:end);
figure
plot(d0,'.-')
hold on
plot(d1,'.-')
grid on; grid minor