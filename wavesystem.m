% z=[0.125,1;0,0;0.1,1]; %test1 data reverse
z=[1,0.125;0,0;1,0.1];
rho=z(1,:); %density
u=z(2,:);%velocity
p=z(3,:);%pressure
y=1.4;
[E]=energy (rho,u,p,y);
pressure=mean(p);
pressure_star=mean(p);
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
u_star=0.5*(fk_right-fk_left)+0.5*(u_right+u_left);

%density
[rho_star_left,rho_star_right]=density_1(p,pressure_star,rho,y);

%shock waves 
[shock_xt_left,shock_xt_right,expansion_xt_left_head,expansion_xt_left_tail,expansion_xt_right_head,expansion_xt_right_tail]=s_e(p,u,pressure_star,y,u_star,rho_star_left,rho_star_right,rho);

p_left=p(1,1);
p_right=p(1,2);
u_left=u(1,1);
u_right=u(1,2);
rho_left=rho(1,1);
rho_right=rho(1,2);
xtL=4;
dxt=0.1;
% if left expansnion 
if isempty(expansion_xt_left_head)==0; 
%set the domain
   xt=[-xtL/2:dxt:u_star];
   n=length(xt);
    for i=1:n
     %after fan area
    if ((xt(i)>=-xtL/2)&&(xt(i)<expansion_xt_left_head));
     x1=(-xtL/2:dxt:expansion_xt_left_head);
%      x1(1,length(x1)) = expansion_xt_left_head;
    end  
     %within fan area
    if ((xt(i)>=expansion_xt_left_head)&&(xt(i)<=expansion_xt_left_tail));
      x2=[expansion_xt_left_head:dxt:expansion_xt_left_tail];   
      [rho_fan_left,u_fan_left,p_fan_left] = fan_left( u,p,y,expansion_xt_left_head,expansion_xt_left_tail,rho,xtL,u_star,dxt);
    elseif(abs(expansion_xt_left_tail-expansion_xt_left_head)<=dxt)
      x2=[];
      rho_fan_left=[];
      u_fan_left=[];
      p_fan_left=[];
    end  
     %before fan area
     if ((xt(i)>expansion_xt_left_tail)&&(xt(i)<=u_star));
      x3=[expansion_xt_left_tail:dxt:u_star];
    end  
    end


  xt_l = [x1,x2,x3];
  rho_l = [rho_left*ones(1,length(x1)),rho_fan_left,rho_star_left*ones(1,length(x3))];
  u_l = [u_left*ones(1,length(x1)),u_fan_left,u_star*ones(1,length(x3))];
  p_l = [p_left*ones(1,length(x1)),p_fan_left,pressure_star*ones(1,length(x3))]; 
%   display('left expansion'); %#ok<DPLAYCHAR>
%----------------------------------------------------------------------------------
%if left shock
elseif isempty(shock_xt_left)==0; 
%set the domain
xt=[-xtL/2:dxt:u_star];
n=length(xt);
 for i=1:n
     if ((xt(i)>=-xtL/2)&&(xt(i)<=shock_xt_left));
         x9=[-xtL/2:dxt:shock_xt_left];
    end  
    %before shock
    if ((xt(i)>=shock_xt_left)&&(xt(i)<u_star));
         x10=[shock_xt_left:dxt:u_star];
    end  
    %after shock
          
 end
  xt_l = [x9,x10];
  rho_l = [rho_left*ones(1,length(x9)),rho_star_left*ones(1,length(x10))];
  u_l = [u_left*ones(1,length(x9)),u_star*ones(1,length(x10))];
  p_l = [p_left*ones(1,length(x9)),pressure_star*ones(1,length(x10))];
%   display('left shock'); %#ok<DPLAYCHAR>
end

  
%---------------------------------------------------------------------------
% if right expansnion 
%set the domain
if isempty(expansion_xt_right_head)==0; 
 xt=[u_star:dxt:xtL/2];
 n=length(xt);
    for i=1:n
      %before fan area
    if ((xt(i)>=u_star)&&(xt(i)<expansion_xt_right_tail));
        x6=[u_star:dxt:expansion_xt_right_tail];
    end  
    %within fan area
    if ((xt(i)>=expansion_xt_right_tail)&&(xt(i)<=expansion_xt_right_head));
        x7=[expansion_xt_right_tail:dxt:expansion_xt_right_head];
        [rho_fan_right,u_fan_right,p_fan_right] = fan_right( u,p,y,expansion_xt_right_head,expansion_xt_right_tail,rho,xtL,u_star,dxt);
    elseif (abs(expansion_xt_right_head-expansion_xt_right_tail)<=dxt)
        x7=[];
        rho_fan_right=[];
        u_fan_right=[];
        p_fan_right=[];
    end     
    %after fan area
     if ((xt(i)>expansion_xt_right_head)&&(xt(i)<=xtL/2));
        x8=[expansion_xt_right_head:dxt:xtL/2];
        x8(1,length(x8)) = xtL/2;
     end
    end
  xt_r = [x6,x7,x8];
  rho_r = [rho_star_right*ones(1,length(x6)),rho_fan_right,rho_right*ones(1,length(x8))];
  u_r = [u_star*ones(1,length(x6)),u_fan_right,u_right*ones(1,length(x8))];
  p_r = [pressure_star*ones(1,length(x6)),p_fan_right,p_right*ones(1,length(x8))]; 
%   display('right expansion'); %#ok<DPLAYCHAR>
%---------------------------------------------------------------------------
%if right shock
elseif isempty(shock_xt_right)==0; 
%set the domain
xt=[u_star:dxt:xtL/2];
n=length(xt);
 for i=1:n
    %before shock
    if ((xt(i)>=u_star)&&(xt(i)<shock_xt_right));
         x4=[u_star:dxt:shock_xt_right];
    end  
    %after shock
    if ((xt(i)>=shock_xt_right)&&(xt(i)<=xtL/2));
        x5=[shock_xt_right:dxt:xtL/2];
        x5(1,length(x5)) = xtL/2;
    end        
 end
  xt_r = [x4,x5];
  rho_r = [rho_star_right*ones(1,length(x4)),rho_right*ones(1,length(x5))];
  u_r = [u_star*ones(1,length(x4)),u_right*ones(1,length(x5))];
  p_r = [pressure_star*ones(1,length(x4)),p_right*ones(1,length(x5))];
%   display('right shock'); %#ok<DPLAYCHAR> 
end
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%information
xt=[xt_l,xt_r];
rho_infor = [rho_l,rho_r];
u_infor = [u_l,u_r];
p_infor = [p_l,p_r];
%---------------------------------------------------------------------------
%find xt=0
xt_interp=xt;
[xt_interp,index]=unique(xt_interp);
xt0_rho=interp1(xt_interp,rho_infor(index),0);
xt0_u=interp1(xt_interp,u_infor(index),0);
xt0_p=interp1(xt_interp,p_infor(index),0);
%---------------------------------------------------------------------------
figure(1)
subplot(1,3,1)
plot(xt,rho_infor);
title('Density');
ylim([0 1.1]);
subplot(1,3,2)
plot(xt,u_infor);
ylim([0 1.1]);
title('velocity');
subplot(1,3,3)
plot(xt,p_infor);
ylim([0 1.1]);
title('Pressure');
