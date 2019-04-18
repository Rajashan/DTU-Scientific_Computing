clear all;
close all;
% deine parmaters
param = [4,1/2];
tspan = [0 10]; % define time interval
h = 0.025; % stepsize
n = (tspan(2)-tspan(1))/h; % number of steps
lambda = 20; % value for lambda
% function of intrst
func = @(t,y) param(1)*t .* y.^param(2) - lambda*(y-(1+t^2)^2); 
y0 = 10; % initial condition
% call our linear Multistep solver
[tout,yout] = LMsolver( func, tspan, n, y0);
% define the solution for the variable z
z = (y0-1)*exp(-lambda*tout);
% Plot results to show the effect of z
figure;
hold on;
plot(tout,yout,'r-*','linewidth',2);
plot(tout,z,'b-','linewidth',2);
legend(strcat('Solution for \lambda=',num2str(lambda),'; y_0=',num2str(y0)),strcat('z(t)=(a-1)exp(-\lambda\cdot t); (\lambda=',num2str(lambda),'; y_0=',num2str(y0),')'),'location','northwest')
xlabel('time t')
set(gca,'fontsize',30)



