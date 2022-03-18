%% Generate Antenna Array
function array = gen_array(n_ants, arrtype, opt, ap, rots)
    
for i = 1:length(ap)
    if arrtype == 1
    % Circular Array
    %
    dst = [];
    ang = 2*pi/n_ants;
    RADIUS = opt.ant_sep/(2 * sin(ang/2));
    for antnum = 1:n_ants
        dst = [dst; 0, RADIUS*sin(ang*(antnum-1)), RADIUS*cos(ang*(antnum-1))];
    end
    array.local_ants{i} = dst;
    array.radii(n_ants) = RADIUS;
    array.RX_ANT_NUM = n_ants;
    %}
    
    elseif arrtype == 2
    % Coprime Array
    %
    dst = [];
    ang_sep = -21.6505/2;
    RADIUS = 0.0770;
    %sep = opt.cent_lambda/2;
    %RADIUS = 4*sep;
    %ang_sep = -rad2deg(atan2(sep/2,sqrt(RADIUS^2-sep^2)));
    angindex = deg2rad([0+ang_sep, 90-ang_sep, 90+ang_sep, 180-ang_sep, 180+ang_sep, 270-ang_sep, 270+ang_sep, 0-ang_sep]);
    for antnum = 1:n_ants
        dst = [dst; 0, RADIUS*sin(angindex(antnum)), RADIUS*cos(angindex(antnum))];
    end
    array.local_ants{i} = dst;
    array.radii(n_ants) = RADIUS;
    array.RX_ANT_NUM = size(dst,1);
    %}
    
    elseif arrtype == 3
    % L Array
    %
    dst = [0, 0, 2*opt.ant_sep;...
            0, 0, 1*opt.ant_sep;...
            0, 0, 0;...
            0, 0, -1*opt.ant_sep;...
            0, 0, -2*opt.ant_sep;...
            0, 0.5*opt.ant_sep, 0;...
            0, -0.5*opt.ant_sep, 0;...
            0, -1.5*opt.ant_sep, 0];

    dst(1:5,3) = dst(1:5,3) - 2*opt.ant_sep;
    dst(6:8,2) = dst(6:8,2) - 1.5*opt.ant_sep;
    
    array.local_ants{i} = dst;
    array.RX_ANT_NUM = 8;
    %}
    
    elseif arrtype == 4
        RX_ANT_NUM = 6;
        dst = [];
        ang = 2*pi/RX_ANT_NUM;
        RADIUS = 2*opt.ant_sep/(2 * sin(ang/2));
        for antnum = 1:RX_ANT_NUM
            dst = [dst; 0, RADIUS*sin(ang*(antnum-1)), RADIUS*cos(ang*(antnum-1))];
        end
        
        RADIUS = sqrt(RADIUS^2-dst(2,3)^2);
        for antnum = 1:RX_ANT_NUM
            dst = [dst; 0, RADIUS*sin(ang/2+ang*(antnum-1)), RADIUS*cos(ang/2+ang*(antnum-1))];
        end
        array.RX_ANT_NUM = 12;
        array.local_ants{i} = dst;
        
    elseif arrtype == 5
    % 3D Axes Array
    %
    if n_ants < 2
        error('Too Few Antennas')
    end
    %ant_sep = 0.736*opt.ant_sep;
    ant_sep = 3/8 * opt.cent_lambda;
    %ant_sep = opt.ant_sep;
    dst = [];
    dst(1,:) = [0,0,0];
    
    for iant = 1:n_ants-1
        dst = [dst; -ant_sep*iant,0,0;...
                    0,-ant_sep*iant,0;...
                    0,0,-ant_sep*iant];
    end   
    dst = dst-mean(dst);
    array.local_ants{i} = dst;
    array.RX_ANT_NUM = (n_ants-1)*3+1;
    %}
    
    elseif arrtype == 6
    % 3D Zig-Zag Axes Array
    %
    if n_ants < 2
        error('Too Few Antennas')
    end
    ang = deg2rad(45);%0.61;
    ant_sep = opt.cent_lambda*3/8 *sqrt(2);
    dst = [];
    dst_x = [0,0,0];dst_y = [0,0,0];dst_z = [0,0,0];
    
    for iant = 1:n_ants-1
        dst_x = [dst_x; -1*iant*ant_sep*cos(ang),0,(-0.5-0.5*(-1)^(iant+1))*ant_sep*sin(ang)];
        
        dst_y = [dst_y; (-0.5-0.5*(-1)^(iant+1))*ant_sep*sin(ang),-1*iant*ant_sep*cos(ang),0];
        
        dst_z = [dst_z; 0,(-0.5-0.5*(-1)^(iant+1))*ant_sep*sin(ang),-1*iant*ant_sep*cos(ang)];
    end   
    dst = [0,0,0;dst_x(2:end,:);dst_y(2:end,:);dst_z(2:end,:)];
    dst = dst-mean(dst);
    array.local_ants{i} = dst;
    array.RX_ANT_NUM = (n_ants-1)*3+1;
    
    elseif arrtype == 7
    % Circular Array with Perpendicular Line through center
    %
    dst = [];
    ang = 2*pi/n_ants;
    RADIUS = opt.ant_sep/(2 * sin(ang/2));
    for antnum = 1:n_ants
        dst = [dst; 0, RADIUS*sin(ang*(antnum-1)), RADIUS*cos(ang*(antnum-1))];
    end
    
    dst = [dst; 0,0,0; opt.ant_sep,0,0; -opt.ant_sep,0,0];
    array.local_ants{i} = dst;
    array.radii(n_ants) = RADIUS;
    array.RX_ANT_NUM = n_ants+3;
    
    elseif arrtype == 8
    % Stacked Circular
    %
    dst = [];
    ang = 2*pi/n_ants;
    RADIUS = opt.ant_sep/(2 * sin(ang/2));
    for antnum = 1:n_ants
        dst = [dst; 0, RADIUS*sin(ang*(antnum-1)), RADIUS*cos(ang*(antnum-1))];
    end
    % Method 1: lambda/2 spacing between two layers
    %h = sqrt(opt.ant_sep.^2 - (2*RADIUS*sin(ang/4)).^2);
    
    % Method 2: Peak Resolution zig-zag
    %res = 1.22*opt.cent_lambda/(2*RADIUS); h = abs(2*RADIUS*sin(ang/4) * tan(res));
    
    % Method 3: "Hexagonal" with lambda/2 spacing
    RADIUS = sqrt(RADIUS.^2 - (opt.ant_sep/2).^2);
    h = sqrt(opt.ant_sep.^2 - (opt.ant_sep/2).^2);
    
    for antnum = 1:n_ants
        dst = [dst; h, RADIUS*sin(ang/2+ang*(antnum-1)), RADIUS*cos(ang/2+ang*(antnum-1))];
    end
    
    %for antnum = 1:n_ants
    %    dst = [dst; 2*h, RADIUS*sin(ang*(antnum-1)), RADIUS*cos(ang*(antnum-1))];
    %end
    
    disp(['Peak Size: ', num2str(rad2deg(1.22*opt.cent_lambda/(2*RADIUS)))])
    disp(['d: ', num2str(norm(dst(1,:)-dst(n_ants+1,:)))])
    disp(['ZZ: ', num2str(rad2deg(asin(h/norm(dst(1,:)-dst(n_ants+1,:)))))]);
    array.local_ants{i} = dst;
    array.radii(n_ants) = RADIUS;
    array.RX_ANT_NUM = size(dst,1);
    
    elseif arrtype == 9
        % Test Array
        load('array_optimization/arrays/test_ant_posdel.mat')
        dst = x(:,1:end).';
        dst = dst-mean(dst);
        array.local_ants{i} = dst;
        array.RX_ANT_NUM = size(dst,1);
        
    elseif arrtype == 10
        % Coprime 4,5 Array
        sep = opt.cent_lambda/2;
        x = unique([4*sep*[0:4], 5*sep*[0:3]]);
        r = 5*sep/sqrt(2);
        cur_d = 0;
        for ii=1:8
            dst(ii,:) = [0,r*sin(x(ii)/r),r*cos(x(ii)/r)];
        end
        dst = [dst; dst+[opt.cent_lambda*3/8,0,0]]; % Stacked Array
        %x = unique([2*sep*[0:2], 3*sep*[0:1]]);
        %tmp = zeros(4,3); tmp(:,1) = x;
        %dst = [dst; tmp];
        array.local_ants{i} = dst;
        array.RX_ANT_NUM = size(dst,1);
    elseif arrtype == 11
        %load('test_ant_pos_1d_opt.mat');
        sep = opt.cent_lambda*4/8;
        %x0 = unique([3*sep*[0:4], 4*sep*[0:3]]);
        %x0 = unique([3*sep*[0:3], 4*sep*[1:5]]);
        x0 = sep*[0:12];
        x = [x0;zeros(1,length(x0))];
        dst = x(1,:).';
        dst = dst-mean(dst);
        dst = [zeros(size(dst,1),1),zeros(size(dst,1),1),dst];
        
        array.local_ants{i} = dst;
        %array.radii(n_ants) = RADIUS;
        array.RX_ANT_NUM = size(dst,1);
    end
    
    rot = rots(i);
    %a = [1,0,ap{i}(1); 0,1,ap{i}(2); 0,0,1];
    %b= [cos(rot),-sin(rot),0; sin(rot),cos(rot),0; 0,0,1];
    %c = [1,0,-ap{i}(1); 0,1,-ap{i}(2); 0,0,1];
    %M(i,:,:) = c'*b*a';
    array.M(i,:,:) = [cos(rot), 0, sin(rot); 0, 1, 0; -sin(rot), 0, cos(rot)];
    for iant = 1:array.RX_ANT_NUM
        %dst(iant,:) = [dst(iant,1:2),1]*squeeze(M(i,:,:))+[0,0,-1+dst(iant,3)];
        dst(iant,:) = squeeze(array.M(i,:,:))*dst(iant,:).';
    end
    dst = dst + ap{i};
    array.ants{i} = dst;
end