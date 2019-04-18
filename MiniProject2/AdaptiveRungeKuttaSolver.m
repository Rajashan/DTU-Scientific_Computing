function [T,X,H,tau,R] = AdaptiveRungeKuttaSolver(fun, tspan, x0, h0, abstol, reltol, varargin)
%Change x0 to avoid errors

nx = length(x0);
x0 = reshape(x0,nx,1);

facmin = 0.1;
facmax = 5;
eps = 0.8;

%Butcher Tableau Information for ERK4C
s  = 4;
A = [0 0 0 0 ; 1/2 0 0 0 ; 0 1/2 0 0 ; 0 0 1 0];
b  = [1/6 ; 1/3 ; 1/3 ; 1/6];
c  = [0; 1/2 ; 1/2 ; 1];

% Initial Parameters and First Output
x  = x0;
h = h0;
t  = tspan(1);
tf = tspan(end);

X = x;
T = t;
H = h0;
tau = 0;
R = 1;

%Allocate Memory for Stage variables
Xs  = zeros(nx,s);
Fs  = zeros(nx,s);

Xr  = zeros(nx,s);
Fr  = zeros(nx,s);

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
        hAs  = h*A;                                         
        hbs  = h*b;
        hcs  = h*c;
        
        Ts = t + hcs; %Time Vector
        % Stage 1
        Xs(:,1) = x;
        Fs(:,1) = feval(fun,Ts(1),Xs(:,1),varargin{:});
        % Stage 2,3,4
        for i=2:s
            Xs(:,i) = x + Fs(:,1:i-1)*hAs(i,1:i-1)';
            Fs(:,i) = feval(fun,Ts(i),Xs(:,i),varargin{:});
        end
        xs = x + Fs*hbs;
        
        %%%%%% Refined Guess Next Iterate %%%%%
        hr = h/2;
        %Refined Step-size Matrix-Vector Products
        hAr = hr*A;                                         
        hbr  = hr*b;
        hcr  = hr*c;
        
        Tr = t + hcr; %Refined Time Vector
        % Stage 1
        Xr(:,1) = x;
        Fr(:,1) = feval(fun,Tr(1),Xr(:,1),varargin{:});
        % Stage 2,3,4
        for i=2:s
            Xr(:,i) = x + Fr(:,1:i-1)*hAr(i,1:i-1)';
            Fr(:,i) = feval(fun,Tr(i),Xr(:,i),varargin{:});
        end
        xr = x + Fr*hbr;
        
        %%% Error Computation %%%
        error = xr - xs;
        r = max(abs(error)./max(abstol,abs(xr).*reltol));
        
        Accept = (r <= 1);
        if Accept
                t = t+h;
                x = xs;
                
                R = [R,r];
                tau = [tau,error(end)];
                H = [H,h];
                T = [T,t];
                X = [X,x];
        end
        %Asympototic Step Size Controller with Limited Change
        h = max( [facmin , min( [facmax , (eps/r)^(1/5) ] ) ] ) * h;
    
    end
end
end