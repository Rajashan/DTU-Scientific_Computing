clear all;
close all;
%% Stability region for the Radau5 method
A = [[(88-7*sqrt(6))/360,(296-169*sqrt(6))/1800,(-2+3*sqrt(6))/225];...
        [(296+169*sqrt(6))/1800,(88+7*sqrt(6))/360,(-2-3*sqrt(6))/225];...
        [(16-sqrt(6))/36,(16+sqrt(6))/36,1/9]]; % Butcher Tableau A matrix
b = [(16-sqrt(6))/36;(16+sqrt(6))/36;1/9]; % Butcher Tableau b vector
e = [1;1;1];

alpha = -25:0.105:25; % interval on the real axis
beta = -25:0.05:25; % inteval on the imaginary axis
nreal = length(alpha);
nimag = length(beta);
% Evaluate R(z) in the selected region of the complex plane
absR = zeros(length(beta),length(alpha));
for kreal = 1:nreal
    for kimag = 1:nimag
        z = alpha(kreal) + i*beta(kimag);
        tmp = inv(eye(size(A))-z*A)*e;
        R = 1 + z*b'*tmp;
        absR(kimag,kreal) = abs(R);
    end
end

%% Plot stability region 
figure;
fs = 20;
map = [0.5, 0.5, 0.5; 1, 1, 1];
hold on;
imagesc(alpha,beta,absR,[0 1]);
h = colorbar;
% contourf(alpha,beta,absR,[0,1])
% colormap(map)
plot([0,0],[-25,25],'k-')
plot([-25,25],[0,0],'k-')
axis image;
axis xy;
xlabel('real');
ylabel('imag');
ylabel(h, '|R(z)|')
title('Stability region for the Radau5')
set(gca,'fontsize',fs)

%% Radau5 on the Van der Pol problem

ta = 0; % starting time
tb = 40; % end time
x0 = [0.5;0.5]; % initial conditions
mu = 10;
param = {mu};

h0 = 10^(-3); % Initial step size
abstol = 10^(-3); % absolute toleranse
reltol = 10^(-3); % relative toleranse
eps = 0.8; % epsilon factor
facmin = 0.5; 
facmax = 5;

func = @(t,x,param) [0, 1; -1, param{1}*(1-x(1)^2)]*[x(1);x(2)]; % define function
Jac = @(t,x,param) [[0,1];[-1-2*param{1}*x(1)*x(2),param{1}*(1-x(1)^2)]];
% Call Radau5 sovler
[tout,yout,stats] = Radau5solver(func,Jac,ta,tb,h0,x0,abstol,reltol,eps,facmin,facmax,param);
% Call matlab solvers ode15s
func = @(t,x) [0, 1; -1, param{1}*(1-x(1)^2)]*[x(1);x(2)]; % define function
opts = odeset('RelTol',reltol,'AbsTol',abstol);
[t15s,y15s] = ode15s(func,[ta tb],x0,opts);

%% Complete plot
figure;
subplot(2,2,1)
hold on;
plot(tout,yout(:,1),'ro','linewidth',2); % Radau5 approximated result
plot(t15s,y15s(:,1),'k-','linewidth',3); % ode15s approximated result
title(strcat('Van der Pol (mu=',num2str(mu),')'))
legend(strcat('Radau5 approximation (',num2str(length(tout)),' steps)'),...
    strcat('ode15s approximation (',num2str(length(t15s)),' steps)'),...
    'Location','northwest')
ylabel('x_1')
xlabel('Time')
xlim([ta tb])
ylim([1.1*min([min(yout(:,1)), min(y15s(:,1))]), 2.5*max([max(yout(:,1)), max(y15s(:,1))])])
set(gca,'fontsize',20)
subplot(2,2,3)
hold on;
plot(tout,yout(:,2),'ro','linewidth',2); % Radau5 approximated result
plot(t15s,y15s(:,2),'k-','linewidth',3); % ode15s approximated result
legend(strcat('Radau5 approximation (',num2str(length(tout)),' steps)'),...
    strcat('ode15s approximation (',num2str(length(t15s)),' steps)'),...
    'Location','northwest')
ylabel('x_2')
xlabel('Time')
xlim([ta tb])
ylim([1.1*min([min(yout(:,2)), min(y15s(:,2))]), 2.5*max([max(yout(:,2)), max(y15s(:,2))])])
set(gca,'fontsize',20)
subplot(2,2,[2,4])
hold on;
plot(yout(:,1),yout(:,2),'ro','linewidth',2); % Radau5 approximated result
plot(y15s(:,1),y15s(:,2),'k-','linewidth',3); % ode15s approximated result
title(strcat('Van der Pol phase plot (mu=',num2str(mu),')'))
legend(strcat('Radau5 approximation (',num2str(length(tout)),' steps)'),...
    strcat('ode15s approximation (',num2str(length(t15s)),' steps)'),...
    'Location','northwest')
ylabel('x_2')
xlabel('x_1')
set(gca,'fontsize',20)

%% plot convergense
figure;
subplot(2,3,1);
plot(tout,stats.hx,'r-o','linewidth',2)
ylabel('Step size, h');
xlabel('Time');
set(gca,'fontsize',20)
subplot(2,3,2);
hold on;
plot(tout,stats.R,'r-','linewidth',2)
plot(tout,stats.R.*0+eps,'k-','linewidth',1.5)
title('Stats from Radau5 solver on the Van der Pol problem')
ylabel('r');
xlabel('Time');
set(gca,'fontsize',20)
subplot(2,3,3);
plot(tout,stats.tau(:,1),'r-','linewidth',2)
ylabel('Error');
xlabel('Time');
set(gca,'fontsize',20)
subplot(2,3,4);
plot([1:1:length(stats.n)],stats.n-1,'r-','linewidth',2)
ylabel('Unaccepted steps, n');
xlabel('Function steps');
set(gca,'fontsize',20)
subplot(2,3,5);
hold on;
plot([1:1:length(stats.n)],stats.k,'k-','linewidth',5)
ylabel('Newton steps, k');
xlabel('Function steps');
set(gca,'fontsize',20)

subplot(2,3,6);
Rejected_steps = sum(stats.n - 1);
Number_of_steps = length(tout);
Function_evaluations = sum(stats.fEval);
Jacobian_evaluations = sum(stats.JEval);
set(gca,'visible','off')
A = table(abstol,reltol,Rejected_steps,Number_of_steps);
B = table(mu,Function_evaluations,Jacobian_evaluations);
% Get the table in string form.
TString1 = evalc('disp(A)');
TString2 = evalc('disp(B)');
% Use TeX Markup for bold formatting and underscores.
TString1 = strrep(TString1,'<strong>','\bf');
TString1 = strrep(TString1,'</strong>','\rm');
TString1 = strrep(TString1,'_','\_');
TString2 = strrep(TString2,'<strong>','\bf');
TString2 = strrep(TString2,'</strong>','\rm');
TString2 = strrep(TString2,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
annotation(gcf,'Textbox','String',TString1,'Interpreter','Tex',...
    'FontName',FixedWidth,'LineStyle','none','Fontsize',15,'Units','Normalized','Position',[0.62 -0.1 0.7 0.5]);
annotation(gcf,'Textbox','String',TString2,'Interpreter','Tex',...
    'FontName',FixedWidth,'LineStyle','none','Fontsize',15,'Units','Normalized','Position',[0.62 -0.22 0.7 0.5]);

