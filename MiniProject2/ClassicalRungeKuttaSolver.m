function [Tout,Xout] = ClassicalRungeKuttaSolver(fun,tspan,x0,h,varargin)
%Change x0 to avoid errors

x0 = reshape(x0,length(x0),1);

%Butcher Tableau Information for ERK4C
s  = 4;
A = [0 0 0 0 ; 1/2 0 0 0 ; 0 1/2 0 0 ; 0 0 1 0];
b  = [1/6 ; 1/3 ; 1/3 ; 1/6];
c  = [0; 1/2 ; 1/2 ; 1];

%Step-size Matrix-Vector Products
hA = h*A;                                         
hb  = h*b;
hc  = h*c;

% Initial Parameters, Step Size, etc.
x  = x0;
t  = tspan(1);
tf = tspan(end);
N = round((tf-t)/h);
nx = length(x0);

%Allocate Memory for Variables
T  = zeros(1,s);
X  = zeros(nx,s);
F  = zeros(nx,s);
Tout = zeros(N+1,1);
Xout = zeros(N+1,nx);

%Initial Values to Output
Tout(1) = t;
Xout(1,:) = x';

%For-Loop for the Stages of each Iteration
for n=1:N
    T = t + hc; %Time Vector
    
    % Stage 1
    X(:,1) = x;
    F(:,1) = feval(fun,T(1),X(:,1),varargin{:});
    
    % Stage 2,3,4
    for i=2:s
        X(:,i) = x + F(:,1:i-1)*hA(i,1:i-1)';
        F(:,i) = feval(fun,T(i),X(:,i),varargin{:});
    end
    
    % Next Iterate
    t = t + h;
    x = x + F*hb;
    
    % Save Iterate to Output
    Tout(n+1) = t;
    Xout(n+1,:) = x';
    
end
end