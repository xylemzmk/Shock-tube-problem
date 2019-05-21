function [rho_fan_left,u_fan_left,p_fan_left] = fan_left( u,p,y,expansion_xt_left_head,expansion_xt_left_tail,rho,xtL,u_star,dxt);
p_left=p(1,1);
u_left=u(1,1);
rho_left=rho(1,1);
a_left=(y*p_left/rho_left)^0.5;

xt=(expansion_xt_left_head:dxt:expansion_xt_left_tail); %1*10 double
xt(1,length(xt))=expansion_xt_left_tail;
n=length(xt);
A=2/(y+1);
B=(y-1)/(y+1);
C=(y-1)/2;
D=(y-1)/(y+1);
E=2*y/(y-1);
rho_fan_left=zeros(1,n);
u_fan_left=zeros(1,n);
p_fan_left=zeros(1,n);
for i=1:n
rho_fan_left(1,i)=rho_left*(A+B*(u_left-xt(1,i))/a_left)^(2/(y-1));
u_fan_left(1,i)=A*(a_left+C*u_left+xt(1,i));
p_fan_left(1,i)=p_left*(A+D*(u_left-xt(1,i))/a_left)^E;
end

