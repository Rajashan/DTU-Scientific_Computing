function [x,stat] = NewtonMethodODE(funjac,tk, xk, h, xinit, tol, maxit, varargin)

%Data storage
stat.converged = false; % converged
stat.nfun      = 0;     % number of function calls
stat.iter      = 0;     % number of iterations
stat.tmp       = 0;     %time stamp
stat.nJac      = 0;     % number of jacobian calls

% Initial values and evaluations
stat.iter = 0;
tk1 = tk + h;
x = xinit;
[f,J] = feval(funjac,tk1,x,varargin{:});
R = x - h*f - xk;
I = eye( size(xk,1) );

while  ( (stat.iter < maxit) && (norm(R,'inf') > tol) )
    stat.iter = stat.iter + 1;      %Update iterations
    H = I - J*h;                    %Define Hessian Approximation
    d = -H\R;                       %Solve for step length
    x = x + d;                      %Calculate new x-value
    
    %Compute function and jacobian at nexfound x-value
    [f,J] = feval(funjac,tk1,x,varargin{:});
    
    %Compute residual function at newfound x-value
    R = x - f*h - xk;
end
end










