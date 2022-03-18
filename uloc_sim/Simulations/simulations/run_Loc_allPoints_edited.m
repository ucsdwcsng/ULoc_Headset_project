function [P_error,fin_err,Loc_estimate] = run_Loc_allPoints_edited(ap_centers,predicted_AoA_all,targetPos)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% This script runs reverse localization on all antennas on a given AP/STA

    % INPUTS
        % cli_pos_for_revloc is file_count x 2 with x and y for each file
        % predicted_AoA is a num_ant_tx x file_count variable with AoA from each tx antenna
        % ap_position is 1 if AP is placed where 201 is placed usually, 2
        % if placed at 202's pos and so on
        % targetPos is the position of the STA (or AP) to be reverse
        % localized
    % OUTPUTS
        % 
        
n_ap = size(predicted_AoA_all,2);
for target_count = 1:size(predicted_AoA_all,1)
    clear x;
    clear y;
    x = 0:1; 
    y = zeros(n_ap,length(x));
%     figure(1),hold on;
    for ap_var = 1:n_ap
        theta_aoa = predicted_AoA_all(target_count,ap_var);  
        
%         plot(targetPos(theta_count,1),targetPos(theta_count,2),'-x');
        
        slope = tan(deg2rad(theta_aoa));
        
        y(ap_var,:) = ap_centers(ap_var,2) + slope * (x - ap_centers(ap_var,1) );
        
        
%         plot(x,y(ap_var,:));
%         scatter( ap_centers(ap_var,1), ap_centers(ap_var,2),'filled');
%         waitforbuttonpress;
    end
%     hold off;

    PA = zeros(n_ap,3); PB  = PA;
    PA(:,1) = x(1,1);
    PA(:,2) = y(:,1);
    PB(:,1) = x(1,end);
    PB(:,2) = y(:,end);
    [P_intersect, ~] = lineIntersect3D(PA, PB);
%     
%     figure(1);
%     plot(targetPos(target_count,1),targetPos(target_count,2),'-x');
%     hold on;
%     plot(x, y);
%     scatter(P_intersect(1),P_intersect(2), 'filled', 'd');
%     hold off;
%     xlabel('X Distance in meters');
%     ylabel('Y Distance in meters');

    Loc_estimate(target_count,:) = P_intersect(1,1:2);
    P_error(target_count,:) = abs((P_intersect(1,1:2) - targetPos(target_count,:)));
    fin_err(target_count) = norm(squeeze(P_error(target_count,2)),squeeze(P_error(target_count,1)));
%     waitforbuttonpress;
end

