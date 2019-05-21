function [xt0_rho,xt0_u,xt0_p] = ifshex(shock_xt_left,shock_xt_right,expansion_xt_left_head,expansion_xt_left_tail,expansion_xt_right_head,expansion_xt_right_tail,u_star,u,rho,p,y,rho_star_left,rho_star_right,pressure_star)
% if left expansnion 
p_left=p(1,1);
p_right=p(1,2);
u_left=u(1,1);
u_right=u(1,2);
rho_left=rho(1,1);
rho_right=rho(1,2);
xtL=6;
dxt=0.01;
if isempty(expansion_xt_left_head)==0; 
%set the domain
   xt=[-xtL/2:dxt:u_star];
   n=length(xt);
    for i=1:n
     %within fan area
    if ((xt(i)>=expansion_xt_left_head)&&(xt(i)<=expansion_xt_left_tail));
      x2=[expansion_xt_left_head:dxt:expansion_xt_left_tail];   
      [rho_fan_left,u_fan_left,p_fan_left] = fan_left( u,p,y,expansion_xt_left_head,expansion_xt_left_tail,rho,xtL,u_star,dxt);
    elseif(abs(expansion_xt_left_tail-expansion_xt_left_head)<=dxt)
      x2=[];
      rho_fan_left=[];
      u_fan_left=[];
      p_fan_left=[];
    end;end
  xt_l = [-xtL/2,expansion_xt_left_head,x2,expansion_xt_left_tail,u_star];
  rho_l = [rho_left,rho_left,rho_fan_left,rho_star_left,rho_star_left];
  u_l = [u_left,u_left,u_fan_left,u_star,u_star];
  p_l = [p_left,p_left,p_fan_left,pressure_star,pressure_star]; 
%   display('left expansion'); %#ok<DPLAYCHAR>
%----------------------------------------------------------------------------------
%if left shock
elseif isempty(shock_xt_left)==0; 
xt_l = [-xtL/2,shock_xt_left,shock_xt_left,u_star];
  rho_l = [rho_left,rho_left,rho_star_left,rho_star_left];
  u_l = [u_left,u_left,u_star,u_star];
  p_l = [p_left,p_left,pressure_star,pressure_star];
%   display('left shock'); %#ok<DPLAYCHAR>
end
%---------------------------------------------------------------------------
% if right expansnion 
%set the domain
if isempty(expansion_xt_right_head)==0; 
 xt=[u_star:dxt:xtL/2];
 n=length(xt);
    for i=1:n
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
    end
  xt_r = [u_star,expansion_xt_right_tail,x7,expansion_xt_right_head,xtL/2];
  rho_r = [rho_star_right,rho_star_right,rho_fan_right,rho_right,rho_right];
  u_r = [u_star,u_star,u_fan_right,u_right,u_right];
  p_r = [pressure_star,pressure_star,p_fan_right,p_right,p_right]; 
%   display('right expansion'); %#ok<DPLAYCHAR>
%---------------------------------------------------------------------------
%if right shock
elseif isempty(shock_xt_right)==0; 
  xt_r = [u_star,shock_xt_right,shock_xt_right,xtL/2];
  rho_r = [rho_star_right,rho_star_right,rho_right,rho_right];
  u_r = [u_star,u_star,u_right,u_right];
  p_r = [pressure_star,pressure_star,p_right,p_right];
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
end

