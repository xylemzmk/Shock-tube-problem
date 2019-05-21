function [fk_left,fk_d_left,fk_right,fk_d_right]=pressure_1(p,pressure,rho,y)
Ak=2./(2.4.*rho);
Bk=p./6;
a=(y.*p./rho).^0.5;
% n=2;
% fk_left=[0,0];
% fk_d_left=[0,0];
% fk_right=[0,0];
% fk_d_right=[0,0];
p_left=p(1,1);
p_right=p(1,2);

%left region
if pressure>p_left %shock
fk_left=(pressure-p(1,1)).*(Ak(1,1)./(pressure+Bk(1,1))).^0.5;
fk_d_left=((Ak(1,1)./(Bk(1,1)+p(1,1))).^0.5).*(1-((pressure-p(1,1))./2.*(Bk(1,1)+p(1,1))));
else   %expansion
fk_left=5*a(1,1).*(((pressure./p(1,1)).^((y-1)/(2*y)))-1);
fk_d_left=((pressure./p(1,1)).^-((y+1)/(2*y)))./(p(1,1).*a(1,1));
end

%right region
if pressure>p_right %shock 
fk_right=(pressure-p(1,2)).*(Ak(1,2)./(pressure+Bk(1,2))).^0.5;
fk_d_right=((Ak(1,2)./(Bk(1,2)+p(1,2))).^0.5).*(1-((pressure-p(1,2))./2.*(Bk(1,2)+p(1,2))));
else   %expansion
fk_right=5*a(1,2).*(((pressure./p(1,2)).^((y-1)/(2*y)))-1);
fk_d_right=((pressure./p(1,2)).^-((y+1)/(2*y)))./(p(1,2).*a(1,2));
end
end
