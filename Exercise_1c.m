clear all;
close all;

k = 2; % Order of derivative

alpha = pi; % Negative starting point of the interval
beta = pi; % End point of the interval
h = pi/15; % Grid step size

xbar = pi/2; % point of interest

x = linspace(-alpha,beta,(alpha+beta)/h); % Grid points for approximation
u = exp(sin(x)); % Original function

f = @(x) exp(sin(x)); % Original function
diff = exp(sin(x)).*cos(x); % First order derivative (analytic)
diff2 = exp(sin(x)).*(cos(x)).^2 + exp(sin(x)).*(-sin(x)); % Second order derivative (analytic)

%% Calculate coefficients c
am = 2; % Number of stencil points in each direction
index = round(3*length(x)/4); % Index of the point of interest in the x-axis
c = fdcoeffF(k,xbar,x(index-am:index+am)); % Calculate coefficients

%% Evaluate derivative

z= zeros(length(x),1); % Creates a vector for the solution to the second derivative
for ii = 3:length(x)-2
    z(ii) = sum(c.*u(ii-2:ii+2)); % % Estimates alle points from the third to the third last
end

%% Return results as text
fprintf(strcat('\n\nTrue value of second derivative in xbar=',num2str(xbar),', d^2 u(xbar)/dx^2=',num2str(diff2(index)),'\n'))
fprintf(strcat('Approximated value of second derivative in xbar=',num2str(xbar),', d^2 u(xbar)/dx^2=',num2str(z(index)),'\n'))
fprintf(strcat('Difference between true solution and approximation tau=',num2str(diff2(index)-z(index)),'\n'))

%% PLOT

figure;
hold on;
fplot(f,'b-','linewidth',2) % Function of interest (analytic)
plot(x,diff2,'r-','linewidth',2) % True second order derivative
plot(x(3:end-2),z(3:end-2),'k--o','linewidth',1.5) % Approximated second order derivative
plot([xbar xbar],[-100,100],'k-.','linewidth',1.2) % plot at vertical line in the point of interest
title('Approximation of the second derivative','fontsize',24)
legend('Original function','Second order diff. (Ana)','Second order diff.(App)','Vertical line at \pi/2','Location','northwest')
xlabel('x')
axis([x(1) x(end) 1.1*min(min(diff2, u)) 1.1*max(max(diff2,u))])
set(gca,'fontsize',24)


