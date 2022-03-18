function DP=compute_multipath_profile2d_fast(h,theta_vals,d_vals,opt)
%     DP=zeros(length(theta_vals),length(d_vals));
%     freq_cent = median(opt.freq);
%     const = 1j*2*pi/3e8;
%     const2 = ij*2*pi*opt.ant_sep*freq_cent/(3e8);
%     h = h.';
%     d_rep = opt.freq'.*(const.*repmat(d_vals,length(opt.freq),1));
%     temp = h*d_rep;
%     theta_rep = (1:size(h,1)).*(const2.*repmat(sin(theta_vals'),1,size(h,1)));
%     DP = theta_rep*temp;
    DP=zeros(length(theta_vals),length(d_vals));
    freqcomp = (1j*2*pi/(3e8)).*opt.freq;
    for i=1:length(theta_vals)
        for j = 1:length(d_vals)
            temp = [1:size(h,2)].*(opt.ant_sep*sin(theta_vals(i)))+d_vals(j);
            coeffecients = transpose(h).*exp(temp'*freqcomp);
            DP(i,j)=sum(coeffecients(:));
        end
    end
end