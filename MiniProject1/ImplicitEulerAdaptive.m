function [T,X,stat] = ImplicitEulerAdaptive(fun,tspan,x0,h0,abstol,reltol,varargin)

%Error Control Parameters
eps = 0.8;
facmin = 0.1;
facmax = 5;

%Root-Finding Parameters
maxit = 100;
tol = 10^(-3);

%Initialization
h = h0;
t0 = tspan(1);
tf = tspan(2);

%Start values
x(:,1) = x0;
t = t0;

%Output
X = x;
T = t;
stat.feval = 0;
stat.iter = 0;
stat.naccept = 0;
stat.nreject = 0;
stat.h = h;
stat.r = 0;

while t < tf
    %Final time step
    if tf < t + h
        h = tf - t; 
    end
    
    %Function Evluation
    f = feval(fun,t,x,varargin{:});
    stat.feval = stat.feval + 1;
    
    AcceptStep = false;
    while ~AcceptStep
     stat.iter = stat.iter + 1;
        
    %Entire Step Start :::
    
    %Initial Guess with Forward Euler
    xe_init = x + h*f;
    %Solve for next iterate
    xe = NewtonMethodODE(fun,t,x,h,xe_init,tol,maxit,varargin{:});
    
    %Entire Step End :::
    
    %Double Step Start :::
    
    %1. Double Step
    hm = 0.5*h;
    %Initial Guess
    xm_init = x + hm*f;
    %Solve for next iterate
    xm = NewtonMethodODE(fun,t,x,hm,xm_init,tol,maxit,varargin{:});
  
    %2. Double Step
    tm = t + hm;
    %Initial Guess
    xf_init = xm + hm * feval(fun,tm,xm,varargin{:});
    stat.feval = stat.feval + 1;
    %Solve for next iterate
    xf = NewtonMethodODE(fun,tm,xm,hm,xf_init,tol,maxit,varargin{:});
    
    %Double Step End :::
    
    %Estimate Error : Entire Step vs Double Step
    e = xf - xe;
    r = max( abs(e) ./  max( abstol , abs(xf).*reltol ) );
    
    %Update if Error O.K otherwise change h and retry
    AcceptStep = ( r <= 1 );
    if AcceptStep
        t = t + h;
        x = xf;
        
        X = [X,x];
        T = [T,t];
        stat.h = [stat.h h];
        stat.r = [stat.r r];
        stat.naccept = stat.naccept + 1;
    else
        stat.nreject = stat.nreject + 1;
    end
    
    %Step Size Controller
    h_update = max( facmin , min( sqrt(eps/r) , facmax ) ) * h;
    h = h_update;
    
    end
    
end
end