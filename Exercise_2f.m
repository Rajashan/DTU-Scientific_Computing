clear all;close all;

% define parameters
param = [4,1/2];
tspan = [0 200]; % time interval

% Define plot color and style
colorStyle = {'r-','b-','m-'};

% Stepsizes, values of lambda and initial conditions
h = [0.025,0.05,0.1];
lambda = [0,1,20]; 
a = [1,2,10];

for jj = 1:3
    figure;
    hold on;
    for ii = 1:3
        b = {};
        for kk = 1:3
            % clear old outputs
            clear yout
            clear tout
            
            % Defines functions of intrest
            func = @(t,y) param(1)*t .* y.^param(2) - lambda(ii)*(y-(1+t^2)^2);
             % number of steps
            n = (tspan(2)-tspan(1))/h(jj); 
            % initial condition
            y0 = a(kk);
            
            % call our Linear Multistep solver
            [tout,yout] = LMsolver( func, tspan, n, y0);
            
            % plot for given value of lambda
            subplot(1,3,ii)
            hold on;
            % plot approximation
            plot(tout,yout,colorStyle{kk},'linewidth',2) 
            
            % define text for legend
            b{length(b)+1} = strcat('Solution for h=',num2str(h(jj)),'; a=',num2str(a(kk)));
        end
        
        % Finishes plot
        title(strcat('Solutions for lambda=',num2str(lambda(ii))));
        legend(b,'fontsize',18,'location','north')
        xlabel('Time t')
        set(gca,'fontsize',20)
        hold off; 
    end
end