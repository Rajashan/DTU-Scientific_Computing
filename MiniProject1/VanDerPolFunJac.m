function [f,J] = VanDerPolFunJac(t,x,mu)

f = zeros(2,1);
f(1,1)=x(2);
f(2,1)=mu*(1-x(1)*x(1))*x(2)-x(1);


J = zeros(2,2);
J(2,1) = -2*mu*x(1)*x(2)-1.0;
J(1,2) = 1.0;
J(2,2) = mu*(1-x(1)*x(1));

