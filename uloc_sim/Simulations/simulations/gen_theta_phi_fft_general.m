function [P,A,dt] = gen_theta_phi_fft_general(H, THETA_VALS, PHI_VALS, opt, ant_pos, plt_profile)

	% H = Complex channel, [N_ant x 1]
	% THETA_VALS, PHI_VALS = search space of theta and phi values, [1 x N_angles], 
	% 						 usually between -90 and 90 deg
	% opt = common signal parameters, like wavelength, frequency etc. 
	% ant_pos = antenna positions in local axes of AP, 3 x N_ant
	% 			for example, for L-antenna ULOC AP, 
	% 			on YZ plane, antennas positions would be 
	% 				[(0, 0, -3l/2), (0, 0, -l), (0, 0, -l/2), (0, 0, 0), 
	% 				 (0 -l/2, 0), (0, -l, 0), (0, -3l/2, 0), (0, -2l, 0)], 
	% 			where l = antenna spacing
    d_vals = 0:0.01:10;
	P = zeros(size(THETA_VALS,2), size(PHI_VALS,2));
    A = zeros(size(THETA_VALS,2), size(PHI_VALS,2));
    
    cir = H*exp(1j*2*pi*d_vals./opt.lambda.');
    [~,taps] = max(abs(cir),[],2);
    taps = floor(median(taps));
    CIR_t = cir(:,taps);
    if plt_profile
        figure(1);subplot(1,2,2)
        plot(d_vals,abs(cir.'))
        title('CIR'); grid on; grid minor
    end
    %CIR = fft((H),size(H,2), 2);
    %[~,taps] = max(abs(CIR),[],2);
    %CIR_t = diag(CIR(:,taps));
    %plot((0:1/256:(1-1/256))*256*3e8/(2*500e6), abs(fftshift(CIR).') )
    Ht = CIR_t;
    
	const = 1j*2*pi./opt.lambda(ceil((end+1)/2));
%     for ii=1:length(THETA_VALS)
%     	for jj=1:length(PHI_VALS)
%     		th = THETA_VALS(ii);
%     		phi = PHI_VALS(jj);
%             wave_vec = repmat(const * [cos(phi)*cos(th); sin(phi); cos(phi)*sin(th)], ...
%                               1, size(ant_pos, 2));              
% 			steering_vec = exp(dot(wave_vec, ant_pos)).';
%     		P(ii, jj) = dot(H, steering_vec);
%     	end
%     end
    tic % Start Timer
    for jj=1:length(PHI_VALS)
        th = THETA_VALS;
        phi = PHI_VALS(jj);
        wave_vec = const.*[cos(phi)*cos(th); sin(phi)*ones(1,length(THETA_VALS)); cos(phi)*sin(th)];            
        steering_vec = exp(wave_vec.' * -ant_pos);%exp(dot(wave_vec, ant_pos)).';
        P(:, jj) = steering_vec*Ht;
        %a = abs(steering_vec.*repmat(H.',length(THETA_VALS),1) + 1);
        %[~,A(:, jj)] = min(a,[],2);
        %A(:,jj,:) = a;
    end
    dt = toc; % End Timer
end