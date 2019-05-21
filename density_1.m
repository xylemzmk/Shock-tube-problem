function [rho_star_left,rho_star_right]=density_1(p,pressure_star,rho,y)
p_left=p(1,1);
p_right=p(1,2);

%left region
if pressure_star>p_left %shock  
A=(pressure_star/p_left)+(y-1)/(y+1);
B=(((y-1)/(y+1))*(pressure_star/p_left))+1;
rho_star_left=rho(1,1)*A/B;
else   %expansion
rho_star_left=rho(1,1)*((pressure_star/p_left)^(1/y));
end

%right region
if pressure_star>p_right %shock  
C=(pressure_star/p_right)+(y-1)/(y+1);
D=(((y-1)/(y+1))*(pressure_star/p_right))+1;
rho_star_right=rho(1,2)*C/D;
else   %expansion
rho_star_right=rho(1,2)*((pressure_star/p_right)^(1/y));
end
end

