%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defines the Steering Matrices given the required search space in theta_val
% and d_vals. For the frequencies defined by the cetral frequency 'f' and
% the subcrier spacing 'df' and antenenna separation 'ant_sep'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = get_2dsteering_matrix(theta_vals,phi_vals,f,ant_sep)    
    THETA = zeros(4,length(theta_vals));
    for itheta = 1:length(theta_vals)
        for iphi = 1:length(phi_vals)
            phi_theta = exp(-1j*2*pi*f/3e8*sin(theta_vals(itheta))*cos(phi_vals(iphi))*ant_sep);
            THETA(1,itheta) = phi_theta; % Ant 1
            THETA(2,itheta) = phi_theta.^2; % Ant 2
            THETA(3,itheta) = phi_theta.^3; % Ant 3
            THETA(4,itheta) = phi_theta.^4; % Ant 4
        end
    end
    
    PHI = zeros(4,length(theta_vals));
    for itheta = 1:length(theta_vals)
        for iphi = 1:length(phi_vals)
            phi_theta = exp(-1j*2*pi*f/3e8*sin(phi_vals(iphi))*ant_sep);
            PHI(1,itheta) = 1; % Ant 1
            PHI(2,itheta) = phi_theta; % Ant 2
            PHI(3,itheta) = phi_theta.^2; % Ant 3
            PHI(4,itheta) = phi_theta.^3; % Ant 3
        end
    end
    
    S = kron(THETA,PHI);
end