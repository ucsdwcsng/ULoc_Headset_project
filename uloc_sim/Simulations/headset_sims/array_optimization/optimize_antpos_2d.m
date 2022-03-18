%% Define global constants
clearvars
%close all
pause('on')
theta_step = deg2rad(3);%0.01;
t_vals = (-deg2rad(90):theta_step:deg2rad(90)).'; % aoa values for FFT
p_vals = (-deg2rad(90):theta_step:deg2rad(90)).';
tp(:,1) = repmat(t_vals, length(p_vals), 1);
tp(:,2) = repmat(p_vals, length(t_vals), 1);

% define opt
c = 3e8;
opt.freq = 4.4928e9;%6489.6e6;
opt.lambda = c./opt.freq;
opt.ant_sep = opt.lambda/2;

W = [cos(tp(:,2)).*cos(tp(:,1)), sin(tp(:,2)), cos(tp(:,2)).*sin(tp(:,1))];

x = [0;0;0];



%box_len = [-5*opt.ant_sep:0.01:5*opt.ant_sep];
box_len = [0*opt.ant_sep:(opt.lambda/10):5*opt.ant_sep];
[ssx, ssy, ssz] = meshgrid(zeros(1,length(box_len)), box_len, box_len);
SS = [ssx(:),ssy(:),ssz(:)];
r_lowbound = 2/8.*opt.lambda;
r_highbound = 4/8.*opt.lambda;

for iN = 1:7
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
figure
S = gen_auto_corr_steering(t_vals.', p_vals.', opt, x);

imagesc(S);

figure
scatter3(x(1,:),x(3,:),x(2,:))
xlabel('x')
ylabel('z')
zlabel('y')
axis equal
%save('test_ant_pos2d.mat', 'x', 'S', '-v7.3')

