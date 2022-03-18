function P = compute_spotfi_profile_vectorized(h, theta_vals, phi_vals, opt,S)
%% Specs
% INPUT 
% h: 3 times 30 matrix of csi measurements (with slope across subcarriers
% removed)
% theta_vals: values of time where the profile has to be evaluate
% d_vals: values of distance where the profile has to be evaluated
% opt: options. opt.freq: frequency array where each frequency corresponds to a subcarrier, opt.ant_sep: antenna separation in m, Planned: threshold
% OUTPUT 
% p: is a length(theta_vals) times length(d_vals) array of complex
% values

%% Steps
% Create the Signal matrix
% Find eigen values
% Detect Noise Subspace
% Compute projects on the noise subspace
% disp('Entered vectorized spotfi');

h=h.';
%{
n_sub = size(h,2);
if(size(h,2)~=n_sub)
    fprintf('h does not have the required number subcarriers. Check the code\n');
end
A = zeros(3*n_sub/2, n_sub);
for i=1:n_sub/2
    A(:,i) = [h(1,i-1+(1:n_sub/2)).';h(2,i-1+(1:n_sub/2)).';h(3,i-1+(1:n_sub/2)).'];
end
for i=n_sub/2+1:n_sub
    A(:,i) = [h(2,i-n_sub/2-1+(1:n_sub/2)).';h(3,i-1-n_sub/2+(1:n_sub/2)).';h(4,i-1-n_sub/2+(1:n_sub/2)).'];
end
R = A*A';
%}
A = h;
R = A*A';
[V, D] = eig(R);
eig_vals= diag(D);
idx = find(eig_vals<opt.threshold*max(abs(eig_vals)));

% [eig_vals,ids] = sort(abs(diag(D)),'descend');
% idx = ids((1 + opt.n_comp):end);

A = 1./vecnorm((V(:,idx))'*S,2).^2;

% Reshape the AoA-ToF Profile matrix
P = reshape(A,length(phi_vals),length(theta_vals));
P = P.';

end