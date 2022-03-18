function h=display_model(model, src, dst, rays, blocked_rays)
    hold on;
    h=plot_polygon(model.walls);
    for i=1:length(h)
        set(h(i),'color','red');
        set(h(i),'linewidth',2);
    end

    for i=1:length(model.reflectors)
        plot(model.reflectors{i}(:,1),model.reflectors{i}(:,2),'g','linewidth',2);
    end

    for i=1:length(model.obstacles)
        h=plot_polygon(model.obstacles{i});
        for j=1:length(h)
            set(h(j),'color','blue');
            set(h(j),'linewidth',2);
        end
    end
    rays = {rays{:}};
    blocked_rays = {blocked_rays{:}};
    for i=1:length(rays) % num_rx x num_tx
        r = rays{i};
        bl = blocked_rays{i};
        for cur_idx=1:length(r) % number of rays between a pair of rx and tx
            ray_to_plot = r{1, cur_idx};
            if bl(cur_idx) == 0
                plot(ray_to_plot(:, 1), ray_to_plot(:, 2), '-b')
            end
        end
    end
    
    scatter([src(:, 1); dst(:, 1)], [src(:, 2); dst(:, 2)]);
    hold off
end