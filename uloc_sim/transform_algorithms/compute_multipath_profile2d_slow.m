function DP=compute_multipath_profile2d_slow(h,theta_vals,d_vals,opt)
    DP=zeros(length(theta_vals),length(d_vals));
    opt.lambda = 3e8./opt.freq;
    for i=1:length(theta_vals)
        for j = 1:length(d_vals)
            for k = 1:size(h,1)
                for l = 1:size(h,2)
                    DP(i,j)=DP(i,j)+h(k,l)*exp(1j*2*pi*(l*opt.ant_sep*sin(theta_vals(i))+d_vals(j))/opt.lambda(k));
                end
            end
        end
    end
%     DP=zeros(length(theta_vals),length(d_vals));
%     for i=1:length(theta_vals)
%         for j = 1:length(d_vals)
%             freqcomp = (1j*2*pi/(3e8)).*opt.freq;
%             temp = [1:size(h,2)].*(opt.ant_sep*sin(theta_vals(i)))+d_vals(j);
%             coeffecients = h.*exp(freqcomp'*temp);
%             DP(i,j)=sum(coeffecients(:));
%         end
%     end
end