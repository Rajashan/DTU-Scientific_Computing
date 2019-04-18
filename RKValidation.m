clear all;
close all;
% Define parameters and the function
param = [4,1/2];
tspan = [0 2]; % time interval
n = 101; % number of time steps
func = @(t,y) param(1)*t .* y.^param(2); % function of intrest
y0 = 1; % initial condition
% Call Runge Kutta sovler
[tout,yout] = RKsolver( func, tspan, n, y0);

% Define the solution from part 2,a
g = tout.^4+2.*tout.^2+1;
%% Plot the results
figure;
hold on;
plot(tout,yout,'r-o','linewidth',2); % Approximated result
plot(tout,g,'b-','linewidth',2); % True solition for lambda=0
title('Solution using the Runge Kutta method')
legend('RK approximation',strcat('Exact solution (\lambda=0)'),'Location','north')
axis([tspan(1) tspan(2) 0 max([max(g), max(yout)])*1.1]);
xlabel('time t')
set(gca,'fontsize',20)

%% Run the code for different stepsizes (h)

for ss = 1:5
    param = [4,1/2];
    tspan = [0 10];  % time interval  
    h(ss) = 1/2^ss; % stepsize
    n = (tspan(2)-tspan(1))/h(ss); % number of grid points
    
    func = @(t,y) param(1)*t .* y.^param(2); % function from part 2,a 
    y0 = 1; % initial condition
    % call Runge Kutta solver
    [tout,yout] = RKsolver( func, tspan, n, y0);
    % Define the solution from part 2,a
    g = tout.^4+2.*tout.^2+1;    
    % Calculate the truncation error by subtracting approximation from the
    % true solution
    tau(ss) = abs(g(round(end/2))-yout(round(end/2))');
end

%% PLot the rate of convergense
figure;
hold on;
plot(h,h.^4,'k-','linewidth',2) % Expected rate of convergense (forth order)
plot(h,tau,'r-.o','linewidth',2) % Estimated error
grid on;
xlabel('step size, h');
ylabel('Error');
title('Rate of convergence (RK)');
legend('Expected curve (fourth order)','Estimated truncation error (RK)','location','northwest');
set(gca, 'XScale', 'log','YScale', 'log','fontsize',20);
xlim([min(h)*0.9 max(h)*1.1])
