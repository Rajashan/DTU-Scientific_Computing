function [f,J] = VanderpolDrift(t,x,p)
% Evaluation of function value and jacobian for the drift term
% in the SDE version of the van der pol. The Jacobian is only evaluated
% if more than 1 output is specified. 

mu = p(1);
tmp = mu*(1.0-x(1)*x(1));
f = zeros(2,1);
f(1,1) = x(2);
f(2,1) = tmp*x(2)-x(1);

if nargout > 1
    J = [0 1; -2*mu*x(1)*x(2)-1.0 tmp];
end