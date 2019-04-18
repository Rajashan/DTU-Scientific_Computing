function [Tout,Xout,H,tau,R] = RungeKutta32AdaptiveSolver(fun, tspan, x0, h0, abstol, reltol, varargin)
%Change x0 to avoid errors

nx = length(x0);
x0 = reshape(x0,nx,1);

facmin = 0.1;
facmax = 5;
eps = 0.8;

%Butcher Tableau Information for Dormand-Prince5(4)
s = 3;
A = [0 0 0 ; 1/2 0 0 ; -1 2 0];
b = [1/6 ; 2/3 ; 1/6];
bhat = [1/4 ; 1/2 ; 1/4];
c = [ 0 ; 1/2 ; 1];
d = b-bhat;

% Initial Parameters and First Output
x  = x0;
h = h0;
t  = tspan(1);
tf = tspan(end);

Xout = x;
Tout = t;
H = h0;
tau = 0;
R = 1;

%Allocate Memory for Stage variables
X  = zeros(nx,s);
F  = zeros(nx,s);

%While-Loop for the Stages of each Iteration
while t < tf
    %Final Stage
    if tf < t + h
        h = tf - t;
    end
    
    Accept = false;
    while ~Accept
        
        %%%%% 1. Guess Next Iterate %%%%%
        
        %Step-size Matrix-Vector Products
        hA  = h*A;                                         
        hb  = h*b;
        hc  = h*c;
        hd  = h*d;
        T = t + hc; %Time Vector
        % Stage 1
        X(:,1) = x;
        F(:,1) = feval(fun,T(1),X(:,1),varargin{:});
        % Stage 2,3,4
        for i=2:s
            X(:,i) = x + F(:,1:i-1)*hA(i,1:i-1)';
            F(:,i) = feval(fun,T(i),X(:,i),varargin{:});
        end
        xs = x + F*hb;
        
        %%% Error Computation %%%
        error = F*hd;
        r = max(abs(error)./max(abstol,abs(x).*reltol));
        
        Accept = (r <= 1);
        if Accept
                t = t+h;
                x = xs;
                
                R = [R,r];
                tau = [tau,error(end)];
                H = [H,h];
                Tout = [Tout,t];
                Xout = [Xout,x];
        end
        %Asympototic Step Size Controller with Limited Change
        h = max( [facmin , min( [facmax , (eps/r)^(1/6) ] ) ] ) * h;
    
    end
end
end