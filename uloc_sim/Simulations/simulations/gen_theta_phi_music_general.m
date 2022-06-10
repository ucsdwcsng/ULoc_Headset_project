function [P,dt] = gen_theta_phi_music_general(H, THETA_VALS, PHI_VALS, opt, ant_pos, n_sigs, plt_profile)

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
    coprime = false;
    if ~coprime
        P = zeros(size(THETA_VALS,2), size(PHI_VALS,2));
        d_vals = 0:0.01:10;
        cir = H*exp(1j*2*pi*d_vals./opt.lambda.');
        [~,taps] = max(abs(cir),[],2);
        taps = floor(median(taps));
        CIR_t = cir(:,taps);
        if plt_profile
            figure(1);subplot(1,2,2)
            plot(d_vals,abs(cir.'))
            title('CIR'); grid on; grid minor
        end
        Ht = CIR_t;
        const = 1j*2*pi/opt.cent_lambda;
        
        tic % Start Timer
        aH = Ht*Ht';%./size(H,2); % TODO: Time average over n_ant times instances, increase rank of aH
        [v,d] = eig(aH);
        [~,idx] = sort(ones(1,size(d,1))*abs(d));
        v = v(:,idx);
        n_sigs = length(find(diag(d)>0.1*max(diag(d))));
        for jj=1:length(PHI_VALS)
            th = THETA_VALS;
            phi = PHI_VALS(jj);
            wave_vec = const * [cos(phi)*cos(th); sin(phi)*ones(1,length(THETA_VALS)); cos(phi)*sin(th)];            
            steering_vec = exp(wave_vec.' * ant_pos);%exp(dot(wave_vec, ant_pos)).';
            P(:,jj) = 1./( sum(abs(conj(steering_vec) * v(:,[1:size(d,1)-n_sigs])).^2, 2) );
        end
        dt = toc; % End Timer
        
%         for jj=1:length(THETA_VALS)
%             th = THETA_VALS(jj);
%             wave_vec = const * [cos(th); 0; sin(th)];            
%             steering_vec = exp(wave_vec.' * ant_pos).';%exp(dot(wave_vec, ant_pos)).';
%             P(jj) = 1./( steering_vec' * v(:,[1:size(d,1)-n_sigs]) * v(:,[1:size(d,1)-n_sigs])' * steering_vec);
%         end
    else
        P = zeros(size(THETA_VALS,2), size(PHI_VALS,2));
        const = 1j*2*pi/opt.lambda;
        m=3; n=4;
        % Determine Unique Spacings
        ds = [];
        for ii = 1:length(ant_pos)
            ds = [ds, ant_pos-ant_pos(:,ii)];
        end
        % Subsample z
        tol = 23e-4; %Max size of dwm1000
        i_A = 1;
        for ii = 2:length(ds)
            if any( vecnorm(ds(:,i_A) - ds(:,ii),2,1)<tol )
                continue
            else
                i_A = [i_A, ii];
            end
        end
        [dss,ii] = sort(ds(3,i_A),2); % Only 1D
        center = ceil(length(dss)./2);
        i_A = i_A(ii(center-m*n:center+m*n));
        %plot(dss,'.-');
        
        Rxx = H*H'; % TODO: Time average over n_ant times instances, increase rank of aH
        z = Rxx(:); z = z(i_A);
        
        Rss = zeros(m*n+1);
        for ii = 1:m*n+1
            Rss = Rss + z((m*n+2-ii):(2*m*n+2-ii))*z((m*n+2-ii):(2*m*n+2-ii))';
        end
        Rss = Rss./(m*n+1);
        virtual_ant_pos = ds(:,i_A(ceil(length(i_A)/2):end));
        
        [v,d] = eig(Rss);
        [~,idx] = sort(ones(1,size(d,1))*abs(d));
        v = v(:,idx);
        %n_sigs = 3;%length(find(d_sorted>0.1*max(d_sorted)));
        for jj=1:length(PHI_VALS)
            th = THETA_VALS;
            phi = PHI_VALS(jj);
            wave_vec = const * [cos(phi)*cos(th); sin(phi)*ones(1,length(THETA_VALS)); cos(phi)*sin(th)];            
            steering_vec = exp(wave_vec.' * virtual_ant_pos);%exp(dot(wave_vec, ant_pos)).';
            P(:,jj) = 1./( sum(abs(conj(steering_vec) * v(:,[1:size(d,1)-n_sigs])).^2, 2) );
        end
        
    end
    
    
end