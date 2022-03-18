function S = gen_auto_corr_steering(THETA_VALS, PHI_VALS, opt, ant_pos)

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

	%P = zeros(size(THETA_VALS,2)* size(PHI_VALS,2), size(ant_pos,2));
	const = 1j*2*pi/opt.lambda;
    P = [];
    t_sv = [];
    p_sv = [];
    for jj=1:length(PHI_VALS)
        th = THETA_VALS;
        phi = PHI_VALS(jj);
        wave_vec = const * [cos(phi)*cos(th); sin(phi)*ones(1,length(THETA_VALS)); cos(phi)*sin(th)];            
        steering_vec = exp(wave_vec.' * -ant_pos);%exp(dot(wave_vec, ant_pos)).';
        P(1+(jj-1)*size(THETA_VALS,2):jj*size(THETA_VALS,2),:) = steering_vec;
    end

%     for ii=1:length(THETA_VALS)+length(PHI_VALS)
%         
%         if ii > length(THETA_VALS)
%             phi = PHI_VALS((ii-length(THETA_VALS)):end);
%             th = THETA_VALS((ii-length(PHI_VALS)):end);
%             t_t = [(ii-length(PHI_VALS)):73];
%             p_t = [(ii-length(THETA_VALS)):37];
%         elseif ii > length(PHI_VALS)
%             phi = PHI_VALS;
%             th = THETA_VALS((1+ii-length(PHI_VALS)):ii);
%             t_t = [(1+ii-length(PHI_VALS)):ii];
%             p_t = [1:37];
%         else
%             phi = PHI_VALS(1:ii);
%             th = THETA_VALS(1:ii);
%             t_t = [1:ii];
%             p_t = [1:ii];
%         end
%         if mod(ii,2) == 0
%             phi = flip(phi);
%             p_t = flip(p_t);
%         else
%             th = flip(th);
%             t_t = flip(t_t);
%         end
%         
%         for jj = 1:length(phi)
%             wave_vec = const * [cos(phi(jj))*cos(th(jj)); sin(phi(jj)); cos(phi(jj))*sin(th(jj))];            
%             steering_vec = exp(wave_vec.' * -ant_pos);
%             P = [P;steering_vec];
%         end
%         p_sv = [p_sv, p_t];
%         t_sv = [t_sv, t_t];
%     end
    S = abs(P*P')./size(ant_pos,2);
end