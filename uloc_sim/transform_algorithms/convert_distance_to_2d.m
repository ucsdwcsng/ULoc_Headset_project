function P_out=convert_distance_to_2d(P,d_vals,d1,d2,ant_pos)
% Treat d1 and d2 as x,y positions.
% ant_pos: 1 by 2 array for the current antennas 
% Supply the antenna positions in the same order as the channels used for
% spotfi code
if(length(P)~=length(d_vals))
    fprintf('Size does not match. Size(P,1) should be length(theta_vals) \n');
end
if(isrow(d1))
    d1=d1';
end
if(isrow(d2))
    d2=d2';
end
if(isrow(d_vals))
    d_vals=d_vals';
end
P_out=zeros(length(d2),length(d1));
X=repmat(d1',length(d2),1);
Y=repmat(d2,1,length(d1));
%T=asin(X./sqrt(Y.^2+X.^2));

%ap_vec = ap_pos(1,:)-ap_pos(end,:);
%T=-cart2pol(X(:),Y(:))+pi/2;
%for i=1:length(X(:))
%    T(i) = sign(sum([X(i),Y(i)].*ap_vec))*(pi/2-acos(abs(sum([X(i),Y(i)].*ap_vec))/norm([X(i),Y(i)])/norm(ap_vec)));
%end
D=sqrt((X-ant_pos(1)).^2+(Y-ant_pos(2)).^2);
%D=sqrt(X.^2+Y.^2);
%[~,T_IDX]=min(abs(repmat(T(:),1,length(theta_vals))-repmat(theta_vals',length(T(:)),1)),[],2);
[~,D_IDX]=min(abs(repmat(D(:),1,length(d_vals))-repmat(d_vals',length(D(:)),1)),[],2);
%IDX=sub2ind(size(P),T_IDX,D_IDX);
P_out(:)=P(D_IDX);
end