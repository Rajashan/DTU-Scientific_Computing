clear all;
close all;
% define parameters and function of intrest
param = [4,1/2];
tspan = [0 2]; % time interval
n = 101; % number of steps
lambda = 0; % Fun outside lambda [-9 25]
func = @(t,y) param(1)*t .* y.^param(2) - lambda*(y-(1+t^2)^2); %function of intrst
y0 = 1; % initial condition
% call our Linear Multistep solver
[tout,yout] = LMsolver( func, tspan, n, y0);

g = tout.^4+2.*tout.^2+1; % Define true solution from part 2,a
%% Plot approximation and true solution
figure;
hold on;
plot(tout,yout,'r-o','linewidth',2) % Results from approximation (LMM)
plot(tout,g,'b-','linewidth',2) % True solition for lambda=0
title(strcat('Solution for \lambda=',num2str(lambda)))
legend('Approximated result',strcat('Exact solution (\lambda=0)'),'Location','north')
axis([tspan(1) tspan(2) 0 max([max(g), max(yout)])*1.1]);
xlabel('time, t')
set(gca,'fontsize',20)

%% Calculate rate of convergense
for ss = 1:5
    % define parameters
    param = [4,1/2]; 
    tspan = [0 10]; % time interval
    h(ss) = 1/2^ss; % stepsize
    n = (tspan(2)-tspan(1))/h(ss); % number of time steps
    func = @(t,y) param(1)*t .* y.^param(2); % function of intrst
    y0 = 1; % initial condition
    % call our Linear Multistep solver
    [tout,yout] = LMsolver( func, tspan, n, y0);
    g = tout.^4+2.*tout.^2+1; % solution for lambda =0
    % Estimate trunctation error by subfrating the approximation from the
    % true solution
    tau(ss) = abs(g(round(end/2))-yout(round(end/2))');
end

%% Plot rate of convergense
figure;
hold on;
plot(h,h.^5,'b-','linewidth',2) % Expected rate of convergense (fifth order)
plot(h,h.^6,'k-','linewidth',2) % Rate of convergense (sixth order)
plot(h,tau,'r-.o','linewidth',2) % Estimated rate of convergense
grid on;
xlabel('step size, h');
ylabel('Error');
title('Rate of convergence');
legend('Expected curve (fifth order)','Sixth order convergence rate','Estimated truncation error','location','northwest');
set(gca, 'XScale', 'log','YScale', 'log','fontsize',20);
xlim([min(h)*0.9 max(h)*1.1])

