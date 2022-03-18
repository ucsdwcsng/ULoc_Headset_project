%% General formulation, L-shaped antenna array
H = [1+1j, 1+2j, 1+3j, 1+4j, 1+5j, 1+6j, 1+7j, 1+8j].';
% H = [0, 0, 0, 0, 0, 0, 0, 0].';
theta_vals = -pi/2:0.01:pi/2;
phi_vals = -pi/2:0.01:pi/2;
l = 2*33.3e-3; % m
opt.lambda = l;
opt.ant_sep = l/2;
ant_pos = [[0, 0, -3*l/2]; [0, 0, -l]; [0, 0, -l/2]; [0, 0, 0]; ...
           [0 -l/2, 0]; [0, -l, 0]; [0, -3*l/2, 0]; [0, -2*l, 0]].';

P_gen = gen_theta_phi_fft_general(H, theta_vals, phi_vals, opt, ant_pos);
figure(1)
imagesc(abs(P_gen))

%% Older specific formulation, L-shaped antenna array
P = compute_theta_phi_fft(H, theta_vals, phi_vals, opt);
figure(2)
imagesc(abs(P))

%%
function P = compute_theta_phi_fft(H, THETA_VALS, PHI_VALS, opt)

    P = zeros(size(THETA_VALS,2));
    theta_length = length(THETA_VALS);
    phi_length = length(PHI_VALS);
    const2 = 1j*2*pi*opt.ant_sep/opt.lambda;

    for p = 1:phi_length
        % NOTE: below line is changed from 1:4 to -1:-1:-4 to follow
        % convention in paper. 
        theta_rep = const2*([-1:-1:-4].*repmat(sin(THETA_VALS.')*cos(PHI_VALS(p)),1,4));
        
        phi_rep = const2*([-3:0].*repmat(ones(theta_length,1)*sin(PHI_VALS(p)),1,4));
        aoa_rep = [phi_rep, theta_rep];
        theta_temp = (exp(aoa_rep)*H);
        %theta_temp = abs(exp(aoa_rep(:,1:4))*H(1:4));
        P(:,p) = theta_temp;
    end
end