% load dataset3_9.mat
% load gfd9aoa.mat
n_ap =4;
n_cli=size(cli_pos,1);
predicted_aoa=zeros(n_ap,n_cli);
actual_aoa=zeros(n_ap,n_cli);

for i=1:n_ap    
    ap_center=mean(ap{i},1);
    for j=1:n_cli
        actual=cli_pos(j,:);
        actual_aoa(i,j)=cart2pol(actual(1)-ap_center(1),actual(2)-ap_center(2));
    end
end

err = (actual_aoa - est_phase') * (180/pi);

 err(err<-180)=err(err<-180)+360;
 err(err>180)=err(err>180)-360;




err = (actual_aoa - est_phase1) * (180/pi);

err = err * (pi/180);
for i = 1:4

    figure;
    plot(actual_aoa(i,:)*(180/pi),'.');
    hold on;
    plot(est_phase(:,i)*(180/pi),'.');
    hold off;
    hold on;
%     plot(cli_pos(:,1),'k');
%     hold off;
%     hold on;
%     plot(cli_pos(:,2),'g');
%     hold off;
    title(i);
    legend('ground truth angle', 'estimated angle')
end

% 
% 
% line(cli_pos(:,1),err(2,:)','Color','r')
% ax1 = gca; % current axes
% ax1.XColor = 'r';
% ax1.YColor = 'r';
% 
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% 
% axis([ax1 ax2],[-2 12 -180 180])
% plot([1:19440],'.','Parent',ax2, 'Color','k')
% 
% 
% figure;
% 
% 
% line(cli_pos(:,2),err(2,:)','Color','r')
% 
% 
% ax3 = gca; % current axes
% ax3.XColor = 'r';
% ax3.YColor = 'r';
% 
% ax3_pos = ax3.Position; % position of first axes
% ax4 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% 
% axis([ax3 ax4],[-1 5 -180 180])
% plot([1:19440],'.','Parent',ax4, 'Color','k')
