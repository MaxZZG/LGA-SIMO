function [outputArg1,outputArg2] = drawConPts(Node)

% Max
% 画出控制点

figure
plot3(Node(1,:),Node(2,:),Node(3,:),'.','color',[1 0 1]);
hold on
% 打印节点编号
txt = 1;
for i = 1 : size(Node,2)
   text(Node(1,i),Node(2,i),Node(3,i),num2str(txt),'FontWeight','bold');
   txt = txt + 1;
   hold on
end


end

