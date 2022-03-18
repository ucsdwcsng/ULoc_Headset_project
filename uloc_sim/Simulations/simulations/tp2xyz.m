function result = tp2xyz(theta,phi,radius, type)
theta = -theta;
phi = -phi;
% type 1 is standard axes, type 2 is generalized axes
if type == 2
    x = radius*sin(deg2rad(phi))*sin(deg2rad(theta));
    y = radius*sin(deg2rad(phi))*cos(deg2rad(theta));
    z = radius*cos(deg2rad(phi));
    result = [x,y,z];
elseif type == 1
    x = radius*sin(deg2rad(phi));
    y = radius*cos(deg2rad(phi))*sin(deg2rad(theta));
    z = -radius*cos(deg2rad(phi))*cos(deg2rad(theta));
    result = [x,y,z];
end