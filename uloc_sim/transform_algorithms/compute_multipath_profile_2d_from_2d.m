function [aoa,tof] = compute_multipath_profile_2d_from_2d(cli_pos,ap)
aoa = zeros(4,1);
tof = zeros(4,1);
for j=1:length(ap)
    ap_pos = mean(ap{j});
    ap_vec=ap{j}(1,:)-ap{j}(end,:);
    X=cli_pos(1)-ap_pos(1);
    Y=cli_pos(2)-ap_pos(2);
    aoa(j)=sign(sum([X,Y].*ap_vec))*(pi/2-acos(abs(sum([X,Y].*ap_vec))/norm([X,Y])/norm(ap_vec)));
    tof(j)=sqrt(X^2+Y^2);
end