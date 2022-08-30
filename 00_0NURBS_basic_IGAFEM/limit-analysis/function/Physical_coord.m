function xloc=Physical_coord(b_net,ni,nj,shl)
global p q

% spit out coordinates;
       icount = 0;
       xloc = zeros(1,2);
       for j = 0:q
           for i = 0:p
             icount = icount + 1;
             xloc(1) = xloc(1) + b_net(ni-i,nj-j,1)*shl(icount);
             xloc(2) = xloc(2) + b_net(ni-i,nj-j,2)*shl(icount);
          end
       end
clear icount
return