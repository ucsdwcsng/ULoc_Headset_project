function channels = compute_channels_from_multipath_profile2d_fast(aoa,tof,strength,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given direct path's AOA and ToF provides the channels corresponding to
% that path alone. Uses 2D FFT to do so

%%
%     H(i,j)=H(i,j)+P(k,l)*exp(-1j*2*pi*(j*opt.ant_sep*sin(theta_vals(k))+d_vals(l))/opt.lambda(i));
%     [m,n] = find(P);
%     H=zeros(length(opt.freq),length(ap{1}));
channels = zeros(length(opt.freq),4,length(aoa));% assuming 4 antennas per AP
% channels_rel = zeros(length(aoa)*(length(aoa-1)),length(opt.freq),4);
for i = 1:length(aoa)
    const1 = -1j*2*pi*opt.ant_sep*sin(aoa)/(3e8);
    const2 = -1j*2*pi*tof/(3e8);
    channels(:,:,i) =  strength.*exp( (const1.*(1:4)'*opt.freq) + (const2.*ones(4,1)*opt.freq) ).';
end
channels = permute(channels,[1,3,2]);
end