clear all;
close all;

k = 2; % Order of derivative

for ss = 2:8
    
    h = 1/2^ss; % Grid step size
    hx(ss-1) = h;
    point = pi/2;

    x = linspace(point-2*h,point+2*h,5); % Grid points for approximation
    u = exp(sin(x)); % Original function

    f = @(x) exp(sin(x)); % Original function
    diff = exp(sin(x)).*cos(x); % First order derivative (analytic)
    diff2 = exp(sin(x)).*(cos(x)).^2 + exp(sin(x)).*(-sin(x)); % Second order derivative (analytic)

    c = fdcoeffF(k,point,x);
    z = sum(c.*u); % Approximation at each point in x

    tau(ss-1) = abs(diff2(3)'-z);

end 
%% PLOT errors

figure;
hold on;
plot([1*10^(-6) hx],[1*10^(-6) hx].^4,'k-','linewidth',3);
plot(hx,tau,'r-o','linewidth',2)
grid on;
xlabel('step size, h');
ylabel('Truncation error, \tau');
title('Rate of convergence');
legend('Reference Line','Error','location','northwest');
set(gca, 'XScale', 'log','YScale', 'log','fontsize',30);
xlim([min(hx)*0.9 max(hx)*1.1])


