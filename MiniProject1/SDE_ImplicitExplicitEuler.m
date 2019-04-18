function [X,stat] = SDE_ImplicitExplicitEuler(fun,gun,T,x0,W,varargin)
%Function for solving the stochastic differential equation using the
%Standard Wiener Process. The function uses an Implicit Euler on the
%deterministic part and an Explicit Euler on the stochastic part

%Memory Allocation
X = nan(length(x0),length(T));
N = length(T) - 1;

%Stats
stat.nfun = 0;  %Function evaluations
stat.lerror = []; %local error
stat.gerror = []; %glocal error

%Time step and initial conditions
X(:,1) = x0;
maxit = 100;
tol = 10^(-8);

%For Loop
for i = 1 : N
    dt = T(:,i+1) - T(:,i);
    dW = W(:,i+1) - W(:,i);

    g  = feval(gun,T(i),X(:,i),varargin{:});
    f = feval(fun,T(i),X(:,i),varargin{:});
    
    %%%%%% Implicit Solution %%%%%%
    
    %Implicit Initial Guess using Explicit Euler
    Psi = X(:,i) + g*dW;
    X_init = f*dt + Psi;

    %Solve next iterate using Newton's Methon
    X(:,i+1) = NewtonMethodSDE(fun,T(i),Psi,dt,X_init,tol,maxit,varargin{:});
    
    %%%%%% Implicit Solution %%%%%%
end
end