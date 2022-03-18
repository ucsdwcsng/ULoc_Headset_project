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

sep = opt.lambda/2;
x0 = unique([4*sep*[0:4], 5*sep*[0:3]]);
x0 = [x0;zeros(1,length(x0))];

x1 = opt.lambda.*1.5*[0:3];
% Define search space SS
%box_len = [-0.095:0.005:0.1];
sdx = 0.005;
[ssx, ssy] = meshgrid(0:sdx:0.5342, 0);
SS = [ssx(:),ssy(:)];

%% Kovarik's Method
if false
    n_ant = 4;
    SSi = SS;
    Ax = W*SSi';
    vecs = exp(1j*2.*pi./opt.lambda.*Ax);
    subject = vecs(:,randi(107,1,n_ant));
    figure(1)
    subplot(2,2,1);imagesc(abs(subject*subject'));axis equal;
    Ak = subject;
    I = eye(361);
    for ii = 1:100
        K = (I-Ak*Ak')*pinv(I+Ak*Ak');
        Ak = (I+K)*Ak;
    end
    subplot(2,2,2);imagesc(abs(Ak*Ak'));axis equal;
    % Extract positions
    th = angle(Ak);
    th = th./(2*pi./opt.lambda);
    locs = pinv(W)*th;
    % Find minimum scaling
    d = locs(1,:);
    df = [];
    for ii = 1:n_ant
        df = [df, d-d(ii)];
    end
    df = unique(abs(df)); df = df(df>1e-3);
    scl = opt.lambda*3/8 ./min(df);
    locs = locs.*scl;
    subplot(2,2,3);
    scatter(locs(1,:),locs(2,:))
    axis equal
    Alocs = exp(1j*2*pi./opt.lambda.*W*conj(locs));
    subplot(2,2,4);
    imagesc(abs(Alocs*Alocs'));axis equal
end

%% Optimizer
SSi = SS;

Ax = W*SSi';
vecs = exp(1j*2.*pi./opt.lambda.*Ax);

y = eye(size(vecs,1));
b = y(:);
for ii = 1:size(vecs,2)
    a = vecs(:,ii)*vecs(:,ii)';
    A2(:,:,ii) = a;
    A(:,ii) = a(:);
    disp([num2str(ii),'/',num2str(size(vecs,2))])
end
disp('All Possible AutoCorr Matrices Found')

%% Exhaustive Method
if false
    counter = 0;
    space = size(A2,3);
    N = 4; % Number Antennas Goal
    selection = [1:N];
    done = false;
    best_err = 1e9;
    while(done == false)
        err = norm(abs(sum(A2(:,:,selection),3))./N-eye(size(A2,1)));
        if err < best_err
            best_err = err;
            best_sel = selection;
        end
        
        for ii = N%N:-1:1
            if selection(ii) < (space-N+ii)
                selection(ii) = selection(ii)+1;
                break
            else
                if ii == 1
                    done = true;
                    break
                else
                    for jj = 1:ii
                        if selection(ii-jj) < (space-N+ii-jj) % If not at cap
                            selection((ii-jj+1):ii) = selection(ii-jj)+1+[1:jj];
                            selection(ii-jj) = selection(ii-jj)+1;
                            break
                        end
                    end
                end
            end
        end % End combination creation
        
        counter = counter + 1;
        if counter == 1000
            disp(selection)
            counter = 0;
        end
        
    end
    
end

%% Least Squares Method
if false
    xx = lsqr(A,b);
    disp('Least Squares Done')
    idx = find(abs(xx)<0.005);
    x = SS(idx,:).';
end
%% Greedy Method
if true % TODO: Try std over autocorr
    x = A(:,1);
    x = A2(:,:,1);
    idx = 1;
    for ii = 1:5
        B = A2 + x;
        %B = sum(abs(B));
        [~,idx(ii+1)] = min(std(B));
        x = x + A(:,idx(ii+1));
    end
    figure()
    imagesc(reshape(abs(x),size(a)))
    x = SS(idx,:).';
end
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
X = exp(1j*2.*pi./opt.lambda.*W*x);
X0 = exp(1j*2.*pi./opt.lambda.*W*x0);
AX = X*X';
AX0 = X0*X0';
figure
subplot(2,2,1);imagesc(abs(AX));title('Optimized 8 Antenna')
subplot(2,2,3);imagesc(abs(AX0));title('Literature 8 Antenna')


subplot(1,2,2)
scatter(x(1,:),x(2,:))
hold on
scatter(x0(1,:),x0(2,:),'r*')
legend('Optimized','Literature')
xlabel('x')
ylabel('z')
zlabel('y')
axis equal

%save('test_ant_pos.mat', 'x', 'S', '-v7.3')

