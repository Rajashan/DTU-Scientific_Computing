function [T,X,stat] = ExplicitEuler(fun,tspan,N,x0,varargin)

%Memory Allocation
X = nan(length(x0),N);
T = nan(1,N);

stat.nfun = 0;  %Function evaluations
stat.lerror = []; %local error
stat.gerror = []; %glocal error

%Time step and initial conditions
dt = ( tspan(2)-tspan(1) ) / N;
X(:,1) = x0;
T(1) = tspan(1);

% Explicit Euler Algorithm with stats. 
for n = 1:N
    X(:,n+1) = X(:,n) + dt*feval(fun,T(n),X(:,n),varargin{:});
    T(n+1) = T(n) + dt;
    
    stat.nfun = stat.nfun + 1;
end

end