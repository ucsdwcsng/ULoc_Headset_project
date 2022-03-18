function DP=compute_multipath_profile(h,ant_pos,lambda,theta_vals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given the N_subcarrierXn_ant channel computes the 1D angle profile
% using 1D FFT
%%
    DP=zeros(size(theta_vals));
    if(iscolumn(ant_pos))
        ant_pos=ant_pos';
    end
    if(iscolumn(h))
        h=h.';
    end
    for i=1:length(theta_vals)
        DP(i)=sum(sum(h.*exp(1j*2*pi*ant_pos*sin(theta_vals(i))./lambda')));
    end
end