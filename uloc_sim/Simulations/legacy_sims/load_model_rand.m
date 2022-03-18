function model = load_model_rand(idx)
%% Define the space
    walls = get_rectangle([0,0],10,8);  % Walls are defined just for plotting. They are non-reflecting, non-blocking                                  
    if(idx==1)
        reflectors{1} = rand(1,2).*[10,8]+[0,0;(rand(1)-0.5)*5,(rand(1)-0.5)*4]; % Reflectors are linear for ease
    elseif (idx==2)
        reflectors{1} = rand(1,2).*[10,8]+[0,0;(rand(1)-0.5)*5,(rand(1)-0.5)*4];
        reflectors{2} = rand(1,2).*[10,8]+[0,0;(rand(1)-0.5)*5,(rand(1)-0.5)*4];
    elseif (idx==3)
        reflectors{1} = rand(1,2).*[10,8]+[0,0;(rand(1)-0.5)*5,(rand(1)-0.5)*4];
        reflectors{2} = rand(1,2).*[10,8]+[0,0;(rand(1)-0.5)*5,(rand(1)-0.5)*4];
        reflectors{3} = rand(1,2).*[10,8]+[0,0;(rand(1)-0.5)*5,(rand(1)-0.5)*4];
    elseif (idx==4)
        reflectors{1} = rand(1,2).*[10,8]+[0,0;(rand(1)-0.5)*5,(rand(1)-0.5)*4];
        reflectors{2} = rand(1,2).*[10,8]+[0,0;(rand(1)-0.5)*5,(rand(1)-0.5)*4];
        reflectors{3} = rand(1,2).*[10,8]+[0,0;(rand(1)-0.5)*5,(rand(1)-0.5)*4];
        reflectors{4} = rand(1,2).*[10,8]+[0,0;(rand(1)-0.5)*5,(rand(1)-0.5)*4];
    elseif (idx==0)
%         obstacles = {};
        reflectors = {};
    
    end
        
    model.walls = walls;
    model.reflectors =  reflectors;
%     model.obstacles = obstacles;
    model.obstacles={};
    %model.reflectors={};
    model.lambda = 3e8./(2.4:0.001:2.48)./1e9;
    model.amps = ones(length(model.lambda),1)*100; % Lambda and amps can be arrays
end