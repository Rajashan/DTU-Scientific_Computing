function [f,J] = GeometricBrownianDrift(t,x,p)

lambda = p(1);
f = lambda*x;

if nargout > 1
    J = lambda;
end