function [T,X,stat] = ImplicitEuler(fun,tspan,N,x0,varargin)

%Memory Allocation
X = nan(length(x0),N);
T = nan(1,N);

%Root Finding Parameters
maxit = 100;
tol = 10^-5;

%Function evaluations and stats
stat.iter = 0;

%Time step and initial conditions
h = ( tspan(2)-tspan(1) ) / N;
X(:,1) = x0;
T(1) = tspan(1);

for n = 1:N
    stat.iter = stat.iter + 1;
    
    %Function evaluation
    X_init = X(:,n) + h*feval(fun,T(n),X(:,n),varargin{:});
    
    %Solve next iterate using Newton's Methon
    X(:,n+1) = NewtonMethodODE(fun,T(n),X(:,n),h,X_init,tol,maxit,varargin{:});
    T(n+1) = T(n) + h;
end
end