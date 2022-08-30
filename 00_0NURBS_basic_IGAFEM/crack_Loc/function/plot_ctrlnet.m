function plot_ctrlnet (CP,col)

% X= reshape (CP(:,:,1), 1, []);
% Y= reshape (CP(:,:,2), 1, []);
% Z= reshape (CP(:,:,3), 1, []);
% 
% plot3 (X,Y,Z,'r.','MarkerSize',20)



  for ii = 1:size (CP, 1)
      X= reshape (CP(ii,:,1), 1, []);
      Y= reshape (CP(ii,:,2), 1, []);
      Z= reshape (CP(ii,:,3), 1, []);
      plot3 (X,Y,Z,'k--')
      plot3 (X,Y,Z,col,'MarkerSize',7,'MarkerFaceColor','r')
%       plot3 (X,Y,Z,'rs ','MarkerSize',10)
  end
    
  for jj = 1:size (CP, 2)
      X= reshape (CP(:,jj,1), 1, []);
      Y= reshape (CP(:,jj,2), 1, []);
      Z= reshape (CP(:,jj,3), 1, []);
      plot3 (X,Y,Z,'k--')
      plot3 (X,Y,Z,col,'MarkerSize',7,'MarkerFaceColor','r')
%       plot3 (X,Y,Z,'rs ','MarkerSize',10)
  end