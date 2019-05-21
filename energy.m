function [E]=energy (xt0_rho,xt0_u,xt0_p,y)
e=xt0_p./((y-1)*xt0_rho); %approimations for the behaviour of e
E=xt0_rho.*(0.5.*xt0_u.^2+e);%energy
end