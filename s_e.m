function [shock_xt_left,shock_xt_right,expansion_xt_left_head,expansion_xt_left_tail,expansion_xt_right_head,expansion_xt_right_tail]=s_e(p,u,pressure_star,y,u_star,rho_star_left,rho_star_right,rho);
p_left=p(1,1);
p_right=p(1,2);
u_left=u(1,1);
u_right=u(1,2);
rho_left=rho(1,1);
rho_right=rho(1,2);
a_left=(y*p_left/rho_left)^0.5;
a_right=(y*p_right/rho_right)^0.5;
a_star_l=(y*pressure_star/rho_star_left)^0.5;
a_star_r=(y*pressure_star/rho_star_right)^0.5;
A=(y+1)/(2*y);
B=(y-1)/(2*y);  


%left region
if pressure_star>p_left %shock
shock_xt_left=u_left-a_left*((A*(pressure_star/p_left)+B)^0.5);
expansion_xt_left_head=[];
expansion_xt_left_tail=[];
else %expansion
shock_xt_left=[];
expansion_xt_left_head=u_left-a_left;
expansion_xt_left_tail=u_star-a_star_l;
end

%right region
if pressure_star>p_right %shock
shock_xt_right=u_right+a_right*(A*(pressure_star/p_right)+B)^0.5;
expansion_xt_right_head=[];
expansion_xt_right_tail=[];
else %expansion
shock_xt_right=[];
expansion_xt_right_head=u_right+a_right;
expansion_xt_right_tail=u_star+a_star_r;
end

end
