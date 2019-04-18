function [X,stat] = SDE_ExplicitExplicitEuler(fun,gun,T,x0,W,varargin)
%Function for solving the stochastic differential equation using the
%Standard Wiener Process.

%Memory Allocation
X = nan(length(x0),length(T));
N = length(T) - 1;


%Stats
stat.nfun = 0;  %Function evaluations
stat.lerror = []; %local error
stat.gerror = []; %glocal error

%Time step and initial conditions
X(:,1) = x0;

%For Loop
for i = 1 : N
    dt = T(:,i+1) - T(:,i);
    dW = W(:,i+1) - W(:,i);

    f = feval(fun,T(i),X(:,i),varargin{:});
    g = feval(gun,T(i),X(:,i),varargin{:});
    
    X(:,i+1) = X(:,i) + f*dt+ g.*dW ;
    stat.nfun = stat.nfun + 1;
end


end