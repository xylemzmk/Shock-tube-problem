function [pressure_star,fk_left,fk_d_left,fk_right,fk_d_right] = pressure_star(u,p,rho,y)
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

end

