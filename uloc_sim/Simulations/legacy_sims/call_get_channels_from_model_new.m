function [channels_wo_offset,channels_w_offset] = call_get_channels_from_model_new(offset_model,no_offset_model,pos,ap,offset)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate channels with and without offsets corresponding to the Models
% offset_model and no_offset_model respectively.

n_ant_per_ap = size(ap{1},1);

channels_wo_offset = zeros(length(no_offset_model.lambda),n_ant_per_ap,length(ap));
channels_w_offset = zeros(length(offset_model.lambda),n_ant_per_ap,length(ap));


for j=1:length(ap)
    
    for k=1:n_ant_per_ap
        
        channels_wo_offset(:,k,j)=get_channels_from_model_edit(no_offset_model,pos,ap{j}(k,:),false,offset(j));
        
        channels_w_offset(:,k,j)=get_channels_from_model_edit(offset_model,pos,ap{j}(k,:),false,0);
        
    end
    
    channels_w_offset(:,:,j) = awgn(squeeze(channels_w_offset(:,:,j)),20);
    channels_wo_offset(:,:,j) = awgn(squeeze(channels_wo_offset(:,:,j)),20);
end
