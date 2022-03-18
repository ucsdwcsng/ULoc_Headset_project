function aoa = get_aoa_for_least_tof(channel,d_vals,theta_vals)
%     aoa = zeros(length(theta_vals),1);
    channels = abs(channel);
    [m,n] = size(channels);
    if(m==length(d_vals))
        channels = channels.';
    elseif(n~=length(d_vals))
        error('Check your channel matrix dimensions')
    end
    bw = imregionalmax(db(channels));
    [ao,tof] = find(bw);
    tof(tof==1) = [];
    tof(tof==length(d_vals)) = [];
    ao(ao==1) = [];
    ao(ao==length(theta_vals)) = [];
    if(isempty(tof) || isempty(ao))
         [~,idx] =  max(abs(channel(:)));
         [ao,tof] = ind2sub(size(channels),idx);
    end
    [~,ind] = min(tof);
    aoa = theta_vals(ao(ind));

end
