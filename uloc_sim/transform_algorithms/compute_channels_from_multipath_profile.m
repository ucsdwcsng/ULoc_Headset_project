function channels_new = compute_channels_from_multipath_profile(P,ant_pos,lambda,theta_vals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given the mulipath profile it converts them back to channels using 1D FFT

%%
    channels_new=zeros(length(lambda),length(ant_pos));
    if(iscolumn(ant_pos))
        ant_pos=ant_pos';
    end
    if(iscolumn(P))
        P=P.';
    end
    if(~iscolumn(theta_vals))
        theta_vals=theta_vals.';
    end
    for i=1:length(lambda)
        channels_new(i,:) = exp(-1j*2*pi*sin(theta_vals(P==1))*ant_pos/lambda(i));
    end
end