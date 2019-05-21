function [rho_fan_right,u_fan_right,p_fan_right] = fan_right( u,p,y,expansion_xt_right_head,expansion_xt_right_tail,rho,xtL,u_star,dxt);
p_right=p(1,2);
u_right=u(1,2);
rho_right=rho(1,2);
a_right=(y*p_right/rho_right)^0.5;

xt=[expansion_xt_right_tail:dxt:expansion_xt_right_head]; %1*10 double
xt(1,length(xt))=expansion_xt_right_head;
n=length(xt);
A=2/(y+1);
B=(y-1)/(y+1);
C=(y-1)/2;
D=(y-1)/(y+1);
E=2*y/(y-1);
rho_fan_right=zeros(1,n);
u_fan_right=zeros(1,n);
p_fan_right=zeros(1,n);
for i=1:n
rho_fan_right(1,i)=rho_right*(A-B*(u_right-xt(1,i)/a_right))^(2/(y-1));
u_fan_right(1,i)=A*(-a_right+C*u_right+xt(1,i));
p_fan_right(1,i)=p_right*(A-D*(u_right-xt(1,i))/a_right)^E;
end

