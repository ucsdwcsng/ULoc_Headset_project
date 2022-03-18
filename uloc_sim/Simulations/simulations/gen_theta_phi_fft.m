function P = gen_theta_phi_fft(H, THETA_VALS, PHI_VALS, opt, theta_angs, phi_angs, ant_sep)

P = zeros(size(THETA_VALS,2));
theta_length = length(THETA_VALS);
phi_length = length(PHI_VALS);
freq_cent = median(opt.freq);
const2 = 1j*2*pi*ant_sep/opt.lambda;

for p = 1:phi_length
    %theta_rep = const2*([1:4].*repmat(sin(THETA_VALS.')*cos(PHI_VALS(p)),1,4));
    %phi_rep = const2*([-3:0].*repmat(ones(theta_length,1)*sin(PHI_VALS(p)),1,4));
    %aoa_rep = [phi_rep, theta_rep];
    aoa_rep = const2.*(sin(repmat(THETA_VALS.',1,8)-theta_angs).*cos(PHI_VALS(p)-phi_angs));
    theta_temp = (exp(aoa_rep)*H.');
    %theta_temp = abs(exp(aoa_rep(:,1:4))*H(1:4));
    P(:,p) = theta_temp;
end

%P(t,p) = P(t,p) + H(antnum)*exp(1j*2*pi*ant_sep(antnum)/opt.lambda.*sin(THETA_VALS(t) - theta_angs(iap,antnum)).*cos(PHI_VALS(p) - phi_angs(iap,antnum)));