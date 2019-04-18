function xdot = preyPredator(t,x,p)
a = p(1);
b = p(2);

%xdot(1) = a * (1 - x(2)) * x(1);
%xdot(2) = -b * (1 - x(1)) * x(2);

xdot = [a * (1 - x(2)) * x(1) ; -b * (1 - x(1)) * x(2)];
%xdot = xdot';

end