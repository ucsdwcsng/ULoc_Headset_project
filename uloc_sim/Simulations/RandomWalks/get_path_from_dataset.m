function walk_idx=get_path_from_dataset(labels_discrete, orients, opt)
    rng;
    walk_pts=opt.walk_pts;
    step_size=opt.step_size;
    margins=opt.margins;    
    walk_idx=zeros(walk_pts,1);    
    
    %Choose a random starting point
    npts=size(labels_discrete,1);
    pt_start=randi(npts);
    walk_idx(1)=pt_start;
    
    %Find a walk
    for i=2:walk_pts
        %Find all points with distance less than step_size
        cur_pt=labels_discrete(walk_idx(i-1),:);
        if ~strcmp(orients, 'none')
            cur_orient=orients(walk_idx(i-1),:);
        end
        d_sq=(cur_pt(1)-labels_discrete(:,1)).^2+(cur_pt(2)-labels_discrete(:,2)).^2;
        idx=find(d_sq<step_size^2 & d_sq>1e-5 & ...
                    labels_discrete(:,1)<max(labels_discrete(:,1))-margins & ...
                    labels_discrete(:,1)>min(labels_discrete(:,1))+margins & ...
                    labels_discrete(:,2)<max(labels_discrete(:,2))-margins & ...
                    labels_discrete(:,2)>min(labels_discrete(:,2))+margins);

        % increase step size until non-zero number of points are within
        % reach
        relax_by=2;
        while(isempty(idx))
            idx=find(d_sq<relax_by*step_size^2 & d_sq>1e-5 & ...
                    labels_discrete(:,1)<max(labels_discrete(:,1))-margins & ...
                    labels_discrete(:,1)>min(labels_discrete(:,1))+margins & ...
                    labels_discrete(:,2)<max(labels_discrete(:,2))-margins & ...
                    labels_discrete(:,2)>min(labels_discrete(:,2))+margins);
            % fprintf('Relaxing by %d for %d \n',relax_by,i);
            relax_by=relax_by*2;
        end
        if strcmp(orients, 'none')
            if(i>2)        
                cur_dir=repmat(cur_pt-labels_discrete(walk_idx(i-2),:),length(idx),1);
                cur_dir=cur_dir./repmat(sqrt(sum(cur_dir.^2,2)),1,2);
                next_dir=labels_discrete(idx,:)-repmat(cur_pt,length(idx),1);
                next_dir=next_dir./repmat(sqrt(sum(next_dir.^2,2)),1,2);
                dir=sum(cur_dir.*next_dir,2);
            else
                dir=ones(length(idx),1);
            end
        else
            cur_dir=repmat([cos(cur_orient), sin(cur_orient)],length(idx),1);
            next_dir= [cos(orients(idx,:)), sin(orients(idx,:))];
            dir=sum(cur_dir.*next_dir,2);
        end
       
        prob=exp(-1*(sqrt(d_sq(idx))-step_size).^2/2/step_size/step_size).*exp(-1*(dir-1).^2*4);
        pt_out=pick_with_probability(idx,prob);
        walk_idx(i)=pt_out;
    end
   
   
end
function pt_out=pick_with_probability(pt_in,prob)
    if(iscolumn(prob))
        prob=prob.';
    end
    prob=prob/sum(prob);
    cum_prob=cumsum(prob);
    
    r=rand();
    if(r<cum_prob(1))
        pt_out=pt_in(1);
    else
        idx=r<cum_prob & r>[0,cum_prob(1:end-1)];
        pt_out=pt_in(idx);
    end

end