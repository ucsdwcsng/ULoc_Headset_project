%% Define global constants
clearvars
%close all
pause('on')
theta_step = deg2rad(30);%0.01;
t_vals = (-deg2rad(180):theta_step:deg2rad(180)).'; % aoa values for FFT
p_vals = (-deg2rad(90):theta_step:deg2rad(90)).';
tp(:,1) = repmat(t_vals, length(p_vals), 1);
tp(:,2) = repmat(p_vals, length(t_vals), 1);

% define opt
c = 3e8;
opt.freq = 4.4928e9;%6489.6e6;
opt.lambda = c./opt.freq;
opt.ant_sep = opt.lambda/2;

W = [cos(tp(:,2)).*cos(tp(:,1)), sin(tp(:,2)), cos(tp(:,2)).*sin(tp(:,1))];

%sep = opt.lambda/2;
%x0 = unique([3*sep*[0:3], 4*sep*[0:2]]);
%x0 = [x0;zeros(1,length(x0))];

% Define search space SS
%box_len = [-0.095:0.005:0.1];
sdx = 0.01;
[ssx, ssy, ssz] = meshgrid(0:sdx:0.2, 0:sdx:0.2, 0:sdx:0.2);
SS = [ssx(:),ssy(:),ssz(:)];

%% Optimizer
SSi = SS;

Ax = W*SSi';
vecs = exp(1j*2.*pi./opt.lambda.*Ax);

for ii = 1:size(vecs,2)
    a = vecs(:,ii)*vecs(:,ii)';
    A(:,ii) = a(:);
    disp([num2str(ii),'/',num2str(size(vecs,2))])
end
disp('Done')

%% Least Squares Method
if false
    b = eye(size(vecs,1));
    b = b(:);
    xx = lsqr(A,b);
    disp('Least Squares Done')
    idx = find(abs(xx)>=3.65e-5);
    x = SS(idx,:).';
end
%% Greedy Method
x = A(:,1);
idx = 1;
for ii = 1:7
    B = A + x;
    B = sum(abs(B));
    [~,idx(ii+1)] = min(B);
    x = x + A(:,idx(ii+1));
    disp(ii)
end
x = SS(idx,:).';
% imagesc(imag(reshape(A(:,1),[181,181])))

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

%% Plots
% figure
% S = gen_auto_corr_steering(t_vals.', p_vals.', opt, x);
% imagesc(S);
%X = exp(1j*2.*pi./opt.lambda.*W*x);
%X0 = exp(1j*2.*pi./opt.lambda.*W*x0);
%AX = X*X';
%AX0 = X0*X0';
S = gen_auto_corr_steering(t_vals.',p_vals.',opt,x);
figure
imagesc(abs(S))
%subplot(2,1,2);imagesc(abs(AX0))


figure
scatter(x(3,:),x(2,:))
hold on
%scatter(x0(1,:),x0(2,:),'r*')
%legend('Custom','Coprime Lit')
xlabel('x')
ylabel('z')
zlabel('y')
axis equal

%save('test_ant_pos.mat', 'x', 'S', '-v7.3')

