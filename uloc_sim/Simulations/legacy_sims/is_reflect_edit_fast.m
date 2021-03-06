function [p,t]=is_reflect_edit_fast(src,dst, reflector)

    % Notes: Only consider reflected paths whose incident rays are close to the
    % normal of the reflector or ones which are incident at the steep angle.
    % For ex: a ray incident at 45 deg will not have a reflection.

    if(sum(abs(src-dst))==0)
        fprintf('SRC and DST are the same in is_reflect. The code does bad things for this case. ');
    end

    % Divide the reflector into 1000 points and check if it reflects
    n_pts = 1000;
    src_angle = zeros(n_pts-1,1);
    dst_angle = zeros(n_pts-1,1);
    % generate the points on the reflector
    cur_point=repmat(reflector(1,:),n_pts-1,1).*repmat((1:(n_pts-1)).',1,size(reflector,2))+repmat(reflector(2,:),n_pts-1,1).*repmat(n_pts-(1:(n_pts-1)).',1,size(reflector,2));
    cur_point=cur_point/n_pts;

    % get the slope of the two parts of the reflector as divided by the
    % cur_point
    line_slope_1=local_cart2pol_fast([reflector(1,1)-cur_point(:,1),reflector(1,2)-cur_point(:,2)]);
    line_slope_2=local_cart2pol_fast([reflector(2,1)-cur_point(:,1),reflector(2,2)-cur_point(:,2)]);

    % get the slope of the lines connecting cur_point with src and dst
    % respectively
    src_slope=local_cart2pol_fast([src(1)-cur_point(:,1),src(2)-cur_point(:,2)]);
    dst_slope=local_cart2pol_fast([dst(1)-cur_point(:,1),dst(2)-cur_point(:,2)]);

    % find the min source angle for each cur_point
    idx=abs(src_slope - line_slope_1)<abs(src_slope-line_slope_2);
    src_angle1(idx)=src_slope(idx)-line_slope_1(idx);
    src_angle1(~idx)=src_slope(~idx)-line_slope_2(~idx);

    % find the min dst angle for each cur_point
    idx=abs(dst_slope - line_slope_1)<abs(dst_slope-line_slope_2);
    dst_angle1(idx) = dst_slope(idx) - line_slope_1(idx);
    dst_angle1(~idx) = dst_slope(~idx) - line_slope_2(~idx);


    [~,idx] = min(abs(src_angle1+dst_angle1));
    t=false; p =[];

    % if the ray just skims the relfector, you want the incident and outgoing
    % ray to be fairly straight
    if(abs(src_angle1(idx)+dst_angle1(idx))<pi)
        p =  reflector(1,:)*idx/n_pts + reflector(2,:)*(n_pts-idx)/n_pts;
        t=true;
    end

    [~,idx] = min(abs(src_angle1-dst_angle1));
    % if the ray is incident close to normal, you want the reflected angle to
    % be close to normal as well
    if((abs(src_angle1(idx)-pi/2) < pi/6 || abs(src_angle1(idx)+pi/2) < pi/6)&&(abs(src_angle1(idx) - dst_angle1(idx)) <pi/3))
        p =  reflector(1,:)*idx/n_pts + reflector(2,:)*(n_pts-idx)/n_pts;
        t=true;
    end

    % find reflector's line equation
    if reflector(1, 1) == reflector(2,1)
        if sign(dst(1) - reflector(1,1)) ~= sign(src(1) - reflector(1,1))
            t = false;
        end
    else
        [coeff,~,~]=polyfit(reflector(:,1),reflector(:,2),1);
        m=coeff(1);
        b=coeff(2);
        if(sign(dst(2)-dst(1)*m-b)~=sign(src(2)-src(1)*m-b))
            %Different sides of reflector
            t=false;
        end
    end

end

function t=local_cart2pol_fast(a)
    t=cart2pol(a(:,1),a(:,2));
end

