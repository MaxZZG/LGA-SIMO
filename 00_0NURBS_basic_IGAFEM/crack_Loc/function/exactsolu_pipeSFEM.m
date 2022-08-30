function [ue,stress_exact]= exactsolu_pipeSFEM(x,y,a,b,E,v,P)
% compute polar info
%a=1;b=5;
    if x<10e-6;
        ang=pi/2;
    else      
        ang=atan(y/x);
    end
    r=sqrt(x^2+y^2);
%ur=(a*a*P*r)*((1-v)+b^2*(1+v)/(r*r))/(E*(b^2-a^2));%for stress plane case
ur=((1+v)*(1-2*v)*(a^2*P)*r^2+(1+v)*a^2*b^2*P)/(E*(b^2-a^2)*r);%for strain plane case
% up(i,1)=((1-v)*(a^2*P)*r(i)^2+(1+v)*a^2*b^2*P)/(E*(b^2-a^2)*r(i));
utheta=0;
ux=ur*cos(ang)-utheta*sin(ang);
uy=ur*sin(ang)+utheta*cos(ang);

ue=[ux;uy];
stressp(1,1)=a^2*(r^2-b^2)*P/((b^2-a^2)*r^2);
stressp(1,2)=a^2*(r^2+b^2)*P/((b^2-a^2)*r^2);
stressp(1,3)=0;

stress_exact(1,1)=stressp(1,1)*cos(ang)^2+stressp(1,2)*sin(ang)^2-0;
stress_exact(2,1)=stressp(1,1)*sin(ang)^2+stressp(1,2)*cos(ang)^2+0;
stress_exact(3,1)=(stressp(1,1)-stressp(1,2))*sin(ang)*cos(ang)+0;


