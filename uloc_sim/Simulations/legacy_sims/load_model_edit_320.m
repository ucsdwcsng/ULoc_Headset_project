function model = load_model_edit(idx)
%% Define the space
    walls = get_rectangle([0,0],9,5);  % Walls are defined just for plotting. They are non-reflecting, non-blocking                                  
    if(idx==1)
        reflectors{1} = [9,1;9,2]; % Reflectors are linear for ease
        reflectors{2} = [9,5;9,4.5]; % Reflectors are linear for ease
        reflectors{3} = [2,4;6,5];
        reflectors{4} = [8,1;2,1];
    elseif (idx==2)
        reflectors{1} = [0,0;2,0];
        reflectors{2} = [0,3;2,5];
        reflectors{3} = [5,5;9,5];
        reflectors{4} = [9,0;6,3];
    elseif (idx==3)
        reflectors{1} = [0,5;5,5];
        reflectors{2} = [0,0;5,2.5];
        reflectors{3} = [7,0.5;7,2];
        reflectors{4} = [7,4;9,4];
    elseif (idx==4)
        reflectors{1} = [0,0;3,1.5];
        reflectors{2} = [0,2;0,5];
        reflectors{3} = [6,1;9,1];
        reflectors{4} = [9,5;5,5];
    elseif (idx==5)
        reflectors{1} = [9,0;9,1];
        reflectors{2} = [0,2.5;0,5];
        reflectors{3} = [0,0;2,0];
        reflectors{4} = [8,2;5,3];
        reflectors{5} = [9,5;5,5];
    elseif (idx==6)
        reflectors{1} = [9,0;9,2];
        reflectors{2} = [0,1;0,2];
        reflectors{3} = [7,2;6,2];
        
    
    end     
    model.walls = walls;
    model.reflectors =  reflectors;
    model.obstacles={};
    freq = double(5e9 + 5*153*1e6) + [-117:-1,1:117].*320e6./256;
    model.lambda = 3e8./freq;
    model.amps = ones(length(model.lambda),1)*100; % Lambda and amps can be arrays
    model.Ts = 3e8/(320e6);
end