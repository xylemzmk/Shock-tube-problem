close all
clear all
L=1000; %domain length
dx=1;
y=1.4;
x=[-L/2-dx:dx:L/2+dx];
N_cell=L/dx+2; %number of cells+ghost
N_face=max(size(x)); %number of faces
infor_cell=zeros(6,N_cell);
infor_cell(1,1:N_cell/2)=0.125; %rho_left
infor_cell(1,N_cell/2+1:N_cell)=1; %rho_right
infor_cell(2,1:N_cell/2)=0; %u_left
infor_cell(2,N_cell/2+1:N_cell)=0; %u_right
infor_cell(3,1:N_cell/2)=0.1; %p_left
infor_cell(3,N_cell/2+1:N_cell)=1; %p_right
%define ghost cells 
infor_cell(1,1)=infor_cell(1,2); infor_cell(1,N_cell)=infor_cell(1,N_cell-1);
infor_cell(2,1)=-infor_cell(2,2); infor_cell(2,N_cell)=-infor_cell(2,N_cell-1);
infor_cell(3,1)=infor_cell(3,2); infor_cell(3,N_cell)=infor_cell(3,N_cell-1);
%initialise Flux array
F=zeros(3,N_face);
%conserved variables
for i=1:N_cell
    infor_cell(4,i)=infor_cell(1,i); %U1
    infor_cell(5,i)=infor_cell(2,i)*infor_cell(1,i); %U2
    infor_cell(6,i)=infor_cell(3,i)/(y-1)+0.5*(infor_cell(2,i)*infor_cell(2,i)*infor_cell(1,i)); %U3
end
%---------------------------------------------------------------------------------------------------
t=0;
tend=7500;
I=0;
record_rho=zeros(tend,N_cell);
while I<=tend
    vmax=max(infor_cell(2,:));
    for i=1:N_cell
      [a] = speedofsound(infor_cell(3,i),infor_cell(1,i),y);
      v=a+abs(infor_cell(2,i));
      if(v >= vmax)
          vmax=v;
      end
    end
    dt=0.5*dx/vmax;
%%
for i=2:N_cell-1  
%i-1/2
rho=[infor_cell(1,i-1),infor_cell(1,i)];
u=[infor_cell(2,i-1),infor_cell(2,i)];
p=[infor_cell(3,i-1),infor_cell(3,i)];
xtL=6;
dxt=0.1;
pressure=0.5*(p(1,1)+p(1,2));
pressure_star=pressure;
error=100;
while error>10^-6
u_left=u(1,1);
u_right=u(1,2);
[fk_left,fk_d_left,fk_right,fk_d_right]=pressure_1(p,pressure,rho,y);
pressure=pressure_star-(fk_left+fk_right+u_right-u_left)/(fk_d_left+fk_d_right);
error=abs((pressure_star-pressure)./0.5.*(pressure_star+pressure));
pressure_star=pressure;
end   
pressure_star=real(pressure_star);
pressure =real(pressure);
%gas velocity
u_left=u(1,1);
u_right=u(1,2);
u_star=0.5*(fk_right-fk_left)+0.5*(u_right+u_left);
%density
[rho_star_left,rho_star_right]=density_1(p,pressure_star,rho,y);
%shock waves 
[shock_xt_left,shock_xt_right,expansion_xt_left_head,expansion_xt_left_tail,expansion_xt_right_head,expansion_xt_right_tail]=s_e(p,u,pressure_star,y,u_star,rho_star_left,rho_star_right,rho);
[xt0_rho,xt0_u,xt0_p] = ifshex(shock_xt_left,shock_xt_right,expansion_xt_left_head,expansion_xt_left_tail,expansion_xt_right_head,expansion_xt_right_tail,u_star,u,rho,p,y,rho_star_left,rho_star_right,pressure_star);
%calculate flux
[E]=energy (xt0_rho,xt0_u,xt0_p,y);
F(1,i)=xt0_rho*xt0_u;
F(2,i)=xt0_rho*xt0_u*xt0_u+xt0_p;
F(3,i)=xt0_u*(E+xt0_p);
%%
%i+1/2
rho=[infor_cell(1,i),infor_cell(1,i+1)];
u=[infor_cell(2,i),infor_cell(2,i+1)];
p=[infor_cell(3,i),infor_cell(3,i+1)];
pressure=0.5*(p(1,1)+p(1,2));
pressure_star=pressure;
error=100;
while error>10^-6
u_left=u(1,1);
u_right=u(1,2);
[fk_left,fk_d_left,fk_right,fk_d_right]=pressure_1(p,pressure,rho,y);
pressure=pressure_star-(fk_left+fk_right+u_right-u_left)/(fk_d_left+fk_d_right);
error=abs((pressure_star-pressure)./0.5.*(pressure_star+pressure));
pressure_star=pressure;
end   
pressure_star=real(pressure_star);
pressure =real(pressure);
%gas velocity
u_left=u(1,1);
u_right=u(1,2);
u_star=0.5*(fk_right-fk_left)+0.5*(u_right+u_left);
%density
[rho_star_left,rho_star_right]=density_1(p,pressure_star,rho,y);
%shock waves 
[shock_xt_left,shock_xt_right,expansion_xt_left_head,expansion_xt_left_tail,expansion_xt_right_head,expansion_xt_right_tail]=s_e(p,u,pressure_star,y,u_star,rho_star_left,rho_star_right,rho);
[xt0_rho,xt0_u,xt0_p] = ifshex(shock_xt_left,shock_xt_right,expansion_xt_left_head,expansion_xt_left_tail,expansion_xt_right_head,expansion_xt_right_tail,u_star,u,rho,p,y,rho_star_left,rho_star_right,pressure_star);
%calculate flux
[E]=energy (xt0_rho,xt0_u,xt0_p,y);
F(1,i+1)=xt0_rho*xt0_u;
F(2,i+1)=xt0_rho*xt0_u*xt0_u+xt0_p;
F(3,i+1)=xt0_u*(E+xt0_p);
%update conserve variable
infor_cell(4:6,i)=infor_cell(4:6,i)+dt/dx*(F(1:3,i)-F(1:3,i+1));
end
%%
%update primitive variable
infor_cell(1,1:N_cell)=infor_cell(4,1:N_cell);
infor_cell(2,1:N_cell)=(infor_cell(5,1:N_cell)./infor_cell(4,1:N_cell));
infor_cell(3,1:N_cell)=(y-1)*(infor_cell(6,1:N_cell)-0.5*(infor_cell(5,1:N_cell).*infor_cell(5,1:N_cell)./infor_cell(4,1:N_cell)));
%update ghost
infor_cell(1,1)=infor_cell(1,2); infor_cell(1,N_cell)=infor_cell(1,N_cell-1);
infor_cell(2,1)=-infor_cell(2,2); infor_cell(2,N_cell)=-infor_cell(2,N_cell-1);
infor_cell(3,1)=infor_cell(3,2); infor_cell(3,N_cell)=infor_cell(3,N_cell-1);
%=============================================================================================================================================
t=t+dt; 
I=I+1
record_rho(I,1:N_cell)=infor_cell(1,1:N_cell);
end
contourf(record_rho,'edgecolor','none');
colormap(gray)
mesh(record_rho)
colormap(gray)

%set('LineStyle','none');
