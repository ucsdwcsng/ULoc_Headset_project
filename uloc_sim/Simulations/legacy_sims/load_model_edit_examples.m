function model = load_model_edit_examples(idx)
%% Define the space
    walls = get_rectangle([0,0],8,8);  % Walls are defined just for plotting. They are non-reflecting, non-blocking                                  
    if(idx==1)
        obstacles = {};
        reflectors = {};
    elseif (idx==2)
        obstacles = {};
        reflectors{1} = [4,0;0,4];
    elseif (idx==3)
        obstacles{1} = get_rectangle([1,1],0.8,0.1);
        reflectors{1} = [4,0;0,4];
    end
        
    model.walls = walls;
    model.reflectors =  reflectors;
    model.obstacles = obstacles;
    freq = double(5e9 + 5*153*1e6) + [-117:-1,1:117].*80e6./256;
    model.lambda = 3e8./freq;
%     model.lambda = 3e8./(2.4:0.001:2.48)./1e9;
    model.amps = ones(length(model.lambda),1)*100; % Lambda and amps can be arrays
    model.Ts = 3e8/(80e6);
end