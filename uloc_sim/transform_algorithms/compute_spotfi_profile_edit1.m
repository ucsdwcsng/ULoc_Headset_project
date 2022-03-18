function P = compute_spotfi_profile_edit1(h, theta_vals, d_vals, opt)
%% Specs
% INPUT 
% h: 4 times 201 matrix of csi measurements (with slope across subcarriers
% removed)
% theta_vals: values of time where the profile has to be evaluate
% t_vals: values of time where the profile has to be evaluated
% opt: options. opt.freq: frequency array where each frequency corresponds to a subcarrier, opt.ant_sep: antenna separation in m, Planned: threshold
% OUTPUT 
% p: is a length(theta_vals) times length(t_vals) array of complex
% values

%% Steps
% Create the Signal matrix
% Find eigen values
% Detect Noise Subspace
% Compute projects on the noise subspace
n_sub = length(opt.freq);
h=h.';
if(size(h,2)~=n_sub)
    fprintf('h does not have 201 subcarriers. Check the code\n');
end
A = zeros(3*n_sub/2, n_sub);
for i=1:(n_sub+1)/2
    A(:,i) = [h(1,i-1+(1:n_sub/2)).';h(2,i-1+(1:n_sub/2)).';h(3,i-1+(1:n_sub/2)).'];
end
for i=n_sub/2+1:n_sub
    A(:,i) = [h(2,i-n_sub/2-1+(1:n_sub/2)).';h(3,i-1-n_sub/2+(1:n_sub/2)).';h(4,i-1-n_sub/2+(1:n_sub/2)).'];
end
% for i=2*n_sub/3+1:n_sub
%     A(:,i) = [h(3,i-n_sub/3-1+(1:n_sub/3)).';h(4,i-1-n_sub/3+(1:n_sub/3)).'];
% end
R = A'*A;
% GR = gpuArray(R);
[V, D] = eig(R);
% V = gather(VR);
% D = gather(DR);
eig_vals= diag(D);
idx = find(eig_vals<0.1*max(abs(eig_vals)));
P = zeros(length(theta_vals), length(d_vals));
sv = zeros(n_sub,length(theta_vals)*length(d_vals));
for i=1:length(theta_vals)
    for j=1:length(d_vals)
        sv(:,(i-1)*length(d_vals)+j) = get_steering_vector(theta_vals(i), d_vals(j), n_sub, mean(opt.freq), mean(diff(opt.freq)), opt.ant_sep);
    end
end
P(:) = diag(1./abs(sv'*V(:,idx)*(V(:,idx))'*sv));
end

function s = get_steering_vector(theta, d, n_sub, f, df, ant_sep)
s = zeros(n_sub,1);
omega_t = exp(-1j*2*pi*df*d/3e8);
phi_theta = exp(-1j*2*pi*f/3e8*sin(theta)*ant_sep);
for i=1:n_sub/2
    s(i) = omega_t^(i-1);
    s(i+n_sub/2) = omega_t^(i-1)* phi_theta;
end

end