% function P = compute_music_profile_2d_fast(h, theta_vals, d_vals, opt)
%% Specs
% INPUT 
% h: 3 times 30 matrix of csi measurements (with slope across subcarriers
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

h = rand(234,4)+1i.*rand(234,4);
theta_vals = -pi/2:0.01:pi/2;
d_vals = -5:0.1:20;
opt.freq = 5e9 + 5*153*1e6 + [-117:-1,1:117].*80e6./256;
opt.ant_sep = 0.0259;
opt.threshold = 0.25;


n_sub = length(opt.freq);
if(size(h,1)~=n_sub)
    fprintf('h does not have 234 subcarriers. Check the code\n');
end
[n_sub,n_ant] = size(h);
A = h(:);
R = A*A';
[V, D] = eig(R);
eig_vals= diag(D);
idx = find(eig_vals<opt.threshold*max(abs(eig_vals)));
P_g = zeros(length(theta_vals).* length(d_vals),1,'gpuArray');
sv = get_steering_vector(theta_vals, d_vals, n_sub, n_ant, mean(opt.freq), mean(diff(opt.freq)), opt.ant_sep);
sv_vec = reshape(sv,size(sv,1),size(sv,2)*size(sv,3));
[sv_len,search_space] = size(sv_vec);
sv_vec_g = gpuArray(sv_vec);
V2 = V(:,idx)*(V(:,idx))';

multiply_search_space = @(x,y) (abs(x'*y*x));

P_g = bsxfun(@mutliply_search_space,sv_vec_g(:,1),V2);

for i=1:search_space
    P_g(i) = 1/abs(squeeze(sv_vec_g(:,i))'*V2*squeeze(sv_vec_g(:,i)));
end

P = reshape(gather(P_g),size(sv,2),size(sv,1));
% end

function d = mutliply_search_space(sv_vec,V2)

d = abs(squeeze(sv_vec)'*V2*squeeze(sv_vec));

end


function s = get_steering_vector(theta, d, n_sub, n_ant, f, df, ant_sep)

omega_t = exp(-1j*2*pi*df*d/3e8);
phi_theta = exp(-1j*2*pi*f/3e8*sin(theta.')*ant_sep);
tau_powers = repmat(0:(n_sub-1),1,n_ant);
phi_powers = repmat(0:(n_ant-1),n_sub,1);
phi_powers = phi_powers(:);
phi_component(:,:,1) = bsxfun( @power, phi_theta, phi_powers.' ).';
tau_component(:,1,:) = bsxfun( @power, omega_t, tau_powers.' );
s =  phi_component.*tau_component;

end