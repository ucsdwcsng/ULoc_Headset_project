function [p,t]=is_reflect_edit_3d(src,dst, reflector)
% Assume quadrilateral reflector defined by 3x3 (n_pts x xyz), 2nd pt is midpoint
% of L defined quadrilateral
if(sum(abs(src-dst))==0)
    fprintf('SRC and DST are the same in is_reflect. The code does bad things for this case. ');
end

% Define cur_point as 10100 points in 3d space along reflector
refline1 = reflector(3,:)-reflector(2,:); refline2 = reflector(1,:)-reflector(2,:);
for ilen = 1:100
    idx = (ilen-1)*101+1;
    cur_point(idx:idx+100,:) = reflector(2,:) + refline1*(ilen-1)/100 ...
                                +[0:0.01:1].'.*refline2;
end

% Find Normal vectors to reflector
norm_vec(1,:) = cross(refline1,refline2); norm_vec(1,:) = norm_vec/norm(norm_vec);
norm_vec(2,:) = -1*norm_vec;

% Define which normal vector to compare to based on minimum distance from
% normal vector coming from center of reflector to src
d1 = norm(src-cur_point(ceil(end/2), :)+norm_vec(1,:));
d2 = norm(src-cur_point(ceil(end/2), :)+norm_vec(2,:));
if d1 > d2
    norm_idx(1) = 2;
else
    norm_idx(1) = 1;
end

% Find vector from each point on reflector to src, then find angle between
% this and the normal
src_vecs = src - cur_point;
src_vecs = src_vecs./el_norm(src_vecs,2);
src_angs = acos( src_vecs*norm_vec(norm_idx(1),:).' ./ el_norm(src_vecs,2) ./ 1);

% Repeat for dst
d1 = norm(dst-cur_point(ceil(end/2), :)+norm_vec(1,:));
d2 = norm(dst-cur_point(ceil(end/2), :)+norm_vec(2,:));
if d1 > d2
    norm_idx(2) = 2;
else
    norm_idx(2) = 1;
end

dst_vecs = dst - cur_point;
dst_vecs = dst_vecs./el_norm(dst_vecs,2);
dst_angs = acos( dst_vecs*norm_vec(norm_idx(2),:).' ./ el_norm(dst_vecs,2) ./ 1);

% Check if same normal was used for both sides, else src/dst on opposite
% sides. If on same side, then find point where src_angs = dst_angs.
if norm_idx(1) ~= norm_idx(2)
    t = false;
    p = []; % reflector blocks
else
    idx = find(abs(src_angs-dst_angs) < 0.1);

    % Now find for which src_vecs and dst_vecs chosen, do their vectors when
    % added and normalized, equal the normal vector
    src_vecs = src_vecs(idx,:); dst_vecs = dst_vecs(idx,:);
    new_vecs = (src_vecs+dst_vecs)./el_norm(src_vecs+dst_vecs,2);
    idx2 = find(el_norm(new_vecs+norm_vec(norm_idx(1),:),2) < 0.05);
    
    if isempty(idx2)
        t = false;
        p = [];
    else
        idx = idx(ceil(median(idx2)));

        t = true;
        p = cur_point(idx,:);
    end
end

end

