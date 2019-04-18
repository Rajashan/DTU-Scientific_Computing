function [f,J] = VanDerPolJac(t,x,p)

f = [x(2) ; p(1) * ( 1 - x(1)^2) * x(2) - x(1)];

J = [0 , -2*p(1)*x(2)*x(1)-1 ; 1 , p(1)*(1-x(1)^2) ];

end