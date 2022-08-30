function plotNURBS_surf(U,V,CP)

mcp = length(CP(:,1,1));
ncp = length(CP(1,:,1));

control=zeros(3,mcp,ncp);
control(1,:,:)=CP(:,:,1);
control(2,:,:)=CP(:,:,2);
control(3,:,:)=CP(:,:,4);


srf = nrbmak(control,{U V});
% nrbctrlplot(srf);
%figure('color',[1 1 1])
% hold on
%figure
nrbkntplot(srf);

