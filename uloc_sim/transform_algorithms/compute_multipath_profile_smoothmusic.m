function P = compute_multipath_profile_smoothmusic(h,theta_vals,d_vals,opt) 
    
    N = size(h,1);
    M = size(h,2);
    Ttot =1;
    K = floor(M/2)+1; 
    L = floor(N/2);
    T = floor(Ttot/2)+1;
    
%     ant_pos = (0:M-1).*opt.ant_sep;
%     if(isrow(ant_pos))
%         ant_pos = ant_pos';       
%     end
    lambda = 3e8./opt.freq;
%     if(iscolumn(h));
%         h=h.';
%     enD
    X = formatCSI(h, N, M, Ttot, K, L, T);
    %H = h'*h;
%     H = X*X';
%     [V, D] = eig(H);   
%     thresh=3;
%     nelem = sum(diag(D) > max(diag(D))/thresh);
%     P=zeros(size(theta_vals));
%     for ii=1:length(theta_vals)
%         e = exp(1j*2*pi*ant_pos*cos(theta_vals(ii))./lambda);
%         P(ii) = 1./abs(e'*V(:, 1:end-nelem)*V(:, 1:end-nelem)'*e);
%     end   
    
    %%
    delayConsider = d_vals;
    u_sConsider = (opt.ant_sep*median(opt.freq)/3e8)*sin( theta_vals );
%     deltaConsider = GridStart(3) + (0:GridPts(3)-1)*GridSpacing(3);
%     delayConsider = unique(delayGridValue, 'stable');
%     u_sConsider = unique(u_sGridValue, 'stable');
%     deltaConsider = unique(deltaGridValue, 'stable');
    
    delaySteeringMat = exp(1j*2*pi*(opt.freq).'*delayConsider);
    aoaSteeringMat = exp(1j*2*pi*((0:K-1)')*u_sConsider);
%     aoaSteeringInvMat = exp(-1i*2*pi*(fc/c)*((0:(T-1))')*deltaConsider);
    
    thetaTauMat = kron(aoaSteeringMat, delaySteeringMat);
%     thetaTauDeltaMat = kron(aoaSteeringInvMat, thetaTauMat);
    H = X*X';
    [V,D] = eig(H);
    thresh=3;
    nelem = sum(diag(D) > max(diag(D))/thresh);
    Qn = V(:,1:end-nelem);
    PnA = (Qn*Qn')*thetaTauMat;
    thetaTauMatTrans = thetaTauMat';
    % we want only the diagonal terms of thetaTauDeltaMatTrans*PnA
    music_spectrum = sum(thetaTauMatTrans.*(PnA.'),2);

    music_spectrum = 1./abs(music_spectrum);
end