function model = load_model_edit(idx)
%% Define the space
    reflectors={};
    walls = get_rectangle([0,0],18,8);  % Walls are defined just for plotting. They are non-reflecting, non-blocking                                  
    if(idx==1)
        reflectors{1} = [16,1;16,2]; % Reflectors are linear for ease
        reflectors{2} = [16,5;16,4.5]; % Reflectors are linear for ease
        reflectors{3} = [4,4;10,5];
        reflectors{4} = [14,1;4,1];
    elseif (idx==2)
        reflectors{1} = [0,0;4,0];
        reflectors{2} = [0,3;4,5];
        reflectors{3} = [10,5;16,5];
        reflectors{4} = [16,0;12,3];
    elseif (idx==3)
        reflectors{1} = [0,5;8,5];
        reflectors{2} = [0,0;8,2.5];
        reflectors{3} = [12,0.5;12,2];
        reflectors{4} = [12,4;16,4];
    elseif (idx==4)
        reflectors{1} = [0,0;6,1.5];
        reflectors{2} = [0,2;0,5];
        reflectors{3} = [10,1;16,1];
        reflectors{4} = [16,5;8,5];
    elseif (idx==5)
        reflectors{1} = [16,0;16,1];
        reflectors{2} = [0,2.5;0,5];
        reflectors{3} = [0,0;2,0];
        reflectors{4} = [14,2;8,3];
        reflectors{5} = [14,5;8,5];
    elseif (idx==6)
        reflectors{1} = [16,0;16,2];
        reflectors{2} = [0,1;0,2];
        reflectors{3} = [12,2;10,2];
    elseif (idx==7)
        reflectors = {};
    end     
    model.walls = walls;
    model.reflectors =  reflectors;
    model.obstacles={};
    freq = double(5e9 + 5*153*1e6) + [-117:-1,1:117].*80e6./256;
    model.lambda = 3e8./freq;
%     model.amps = ones(length(model.lambda),1)*100; % Lambda and amps can be arrays
    model.amps = ones(length(model.lambda),1)*1000;
    model.Ts = 3e8/(80e6);
end