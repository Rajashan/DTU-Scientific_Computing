function [W,Tw,dW,S] = MultiWienerProcess(T,N,nW,Ns,seed)
%Multivariate Standard Brownian Motion = Standard Wiener Process

%Control RNG's (randn)
if nargin == 5
rng(seed)
end

dt = T/N;
%Each row in the 2D matrix represents movement along some axis.
%Each column is the random movement
%Each 3D is a new particle/realization
dW = sqrt(dt) * randn(nW,N,Ns);
W = [ zeros(nW,1,Ns) , cumsum(dW,2)];
Tw = 0:dt:T;

%Compute means by the mean of all realizations
u = mean(W,3);
s = std(W,0,3);

%The 3D matrix is reduced to 2D (Rows are dimensions, Columns are
%time-steps)
u_p2s = u + 2*s;
u_m2s = u - 2*s;

S = [u ; u_p2s ; u_m2s ];
end

