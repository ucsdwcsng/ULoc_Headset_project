clear
for i=1:10
    fn=sprintf('dataset_distance_aoa_large_field_%d.mat',i);
    load(fn);
    max_x =max(cur_model.walls(:,1));
    max_y =max(cur_model.walls(:,2));
    d1 = -max_x:1:2*max_x;
    d2=-max_y:1:2*max_y;
    labels_gaussian_2d=get_gaussian_labels(labels((i-1)*10000+1:i*10000,:),2,d1,d2);
    save(fn,'ap','cur_model','features','labels','labels_gaussian_2d','-v7.3');
    disp(i);
end
