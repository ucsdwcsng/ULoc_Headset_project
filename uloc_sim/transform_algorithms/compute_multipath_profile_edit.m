function DP=compute_multipath_profile_edit(h,ant_pos,lambda,theta_vals,d_vals)
    DP=zeros(length(theta_vals),length(lambda));
    if(iscolumn(ant_pos))
        ant_pos=ant_pos';
    end
    if(iscolumn(h))
        h=h.';
    end
    for i=1:length(theta_vals)
        for j = 1:length(lambda)
            DP(i,j)=sum(h(j).*exp(1j*2*pi*ant_pos*sin(theta_vals(i))./lambda(j)));
        end
    end
end