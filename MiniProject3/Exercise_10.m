clear all;
close all;
%% Test equation

ta = 0; % Starting time
tb = 2; % end time
x0 = 1; % Initial condition
lambda = 1; 
param = {lambda};
func = @(t,x,param) lambda*x; % function of intrest (test equation)
Jac = @(t,x,param) param{1};

h0 = 10^(-3); % Initial step size
abstol = 10^(-3); % absolute toleranse
reltol = 10^(-3); % relative toleranse
eps = 0.8; % epsilon factor
facmin = 0.1; 
facmax = 5;

[tout,yout,stats] = ESDIRK23solver(func,Jac,ta,tb,h0,x0,abstol,reltol,eps,facmin,facmax,param);

g = x0*exp(lambda*tout); % Exact solution when ta=0

figure;
hold on;
plot(tout,g,'k-','linewidth',4);
plot(tout,yout,'r--','linewidth',4);
title('ESDIRK23 solution to the test equation');
legend('Exact solution','ESDIRK23 solution','location','north');
xlabel('Time');
xlim([ta tb]);
ylim([0 max([max(yout), max(g)])]);
set(gca,'fontsize',20);

% convergense plots
figure;
subplot(2,3,1);
plot(tout,stats.hx,'r-o','linewidth',2)
ylabel('Step size, h');
xlabel('Time');
set(gca,'fontsize',20)
subplot(2,3,2);
plot(tout,stats.tau(:,1),'r-','linewidth',2)
title('Stats from ESDIRK solver on the test equation')
ylabel('Error');
xlabel('Time');
set(gca,'fontsize',20)
subplot(2,3,3);
hold on;
plot(tout,stats.R,'r-','linewidth',2)
plot(tout,stats.R.*0+eps,'k-','linewidth',1.5)
ylabel('r');
xlabel('Time');
set(gca,'fontsize',20)
subplot(2,3,4);
plot([1:1:length(stats.n)],stats.n-1,'r-','linewidth',2)
ylabel('Unaccepted steps, n');
xlabel('Function steps');
set(gca,'fontsize',20)
subplot(2,3,5);
hold on;
plot([1:1:length(stats.n)],stats.k(:,1),'r-','linewidth',2)
plot([1:1:length(stats.n)],stats.k(:,2),'b-','linewidth',1.5)
ylabel('Newton steps, k');
xlabel('Function steps');
set(gca,'fontsize',20)

%% Van der Pol problem

ta = 0;
tb = 50;
x0 = [0.5;0.5]; % initial conditions
mu = 10;
param = {mu};

h0 = 10^(-3); % Initial step size
abstol = 10^(-5); % absolute toleranse
reltol = 10^(-5); % relative toleranse
eps = 0.8; % epsilon factor
facmin = 0.5; 
facmax = 5;

func = @(t,x,param) [0, 1; -1, param{1}*(1-x(1)^2)]*[x(1);x(2)]; % define function
Jac = @(t,x,param) [[0,1];[-1-2*param{1}*x(1)*x(2),param{1}*(1-x(1)^2)]];
% Call RK3(2) sovler
[tout,yout,stats] = ESDIRK23solver(func,Jac,ta,tb,h0,x0,abstol,reltol,eps,facmin,facmax,param);
% Call matlab solvers ode15s
func = @(t,x) [0, 1; -1, param{1}*(1-x(1)^2)]*[x(1);x(2)]; % define function
opts = odeset('RelTol',reltol,'AbsTol',abstol);
[t15s,y15s] = ode15s(func,[ta tb],x0,opts);

% Complete plot
figure;
subplot(2,2,1)
hold on;
plot(tout,yout(:,1),'ro','linewidth',2); % ESDIRK23 approximated result
plot(t15s,y15s(:,1),'k-','linewidth',2); % ode15s approximated result
title(strcat('Van der Pol (mu=',num2str(mu),')'))
legend(strcat('ESDIRK23 approximation (',num2str(length(tout)),' steps)'),...
    strcat('ode15s approximation (',num2str(length(t15s)),' steps)'),...
    'Location','northwest')
ylabel('x_1')
xlabel('Time')
xlim([ta tb])
ylim([1.1*min([min(yout(:,1)), min(y15s(:,1))]), 2.5*max([max(yout(:,1)), max(y15s(:,1))])])
set(gca,'fontsize',20)
subplot(2,2,3)
hold on;
plot(tout,yout(:,2),'ro','linewidth',2); % ESDIRK23 approximated result
plot(t15s,y15s(:,2),'k-','linewidth',2); % ode15s approximated result
legend(strcat('ESDIRK23 approximation (',num2str(length(tout)),' steps)'),...
    strcat('ode15s approximation (',num2str(length(t15s)),' steps)'),...
    'Location','northwest')
ylabel('x_2')
xlabel('Time')
xlim([ta tb])
ylim([1.1*min([min(yout(:,2)), min(y15s(:,2))]), 2.5*max([max(yout(:,2)), max(y15s(:,2))])])
set(gca,'fontsize',20)
subplot(2,2,[2,4])
hold on;
plot(yout(:,1),yout(:,2),'ro','linewidth',2); % ESDIRK23 approximated result
plot(y15s(:,1),y15s(:,2),'k-','linewidth',3); % ode15s approximated result
title(strcat('Van der Pol phase plot (mu=',num2str(mu),')'))
legend(strcat('ESDIRK23 approximation (',num2str(length(tout)),' steps)'),...
    strcat('ode15s approximation (',num2str(length(t15s)),' steps)'),...
    'Location','northwest')
ylabel('x_2')
xlabel('x_1')
set(gca,'fontsize',20)

%% plot convergense
figure;
subplot(2,3,1);
plot(tout,stats.tau(:,1),'r-','linewidth',2)
ylabel('Error');
xlabel('Time');
set(gca,'fontsize',20)
subplot(2,3,2);
plot(tout,stats.hx,'r-o','linewidth',2)
title('Stats from ESDIRK solver on the Van der Pol problem')
ylabel('Step size, h');
xlabel('Time');
set(gca,'fontsize',20)
subplot(2,3,3);
hold on;
plot(tout,stats.R,'r-','linewidth',2)
plot(tout,stats.R.*0+eps,'k-','linewidth',1.5)
ylabel('r');
xlabel('Time');
set(gca,'fontsize',20)
subplot(2,3,4);
plot([1:1:length(stats.n)],stats.n-1,'r-','linewidth',2)
ylabel('Unaccepted steps, n');
xlabel('Function steps');
set(gca,'fontsize',20)
subplot(2,3,5);
hold on;
plot([1:1:length(stats.n)],stats.k(:,1),'r-','linewidth',2)
plot([1:1:length(stats.n)],stats.k(:,2),'b-','linewidth',1)
legend('Second inner stage, X_2','Third inner stage, X_3','location','north')
ylabel('Newton steps, k');
xlabel('Function steps');
set(gca,'fontsize',20)

subplot(2,3,6);
Rejected_steps = sum(stats.n - 1);
Newton_interations = sum(stats.k(:,1));
Mean_Newton_interations = mean(stats.k(:,1));
Number_of_steps = length(tout);
Function_evaluations = sum(stats.fEval);
Jacobian_evaluations = sum(stats.JEval);
set(gca,'visible','off')
A = table(abstol,reltol,Rejected_steps,Number_of_steps);
B = table(Newton_interations,Mean_Newton_interations);
C = table(Function_evaluations,Jacobian_evaluations);
% Get the table in string form.
TString1 = evalc('disp(A)');
TString2 = evalc('disp(B)');
TString3 = evalc('disp(C)');
% Use TeX Markup for bold formatting and underscores.
TString1 = strrep(TString1,'<strong>','\bf');
TString1 = strrep(TString1,'</strong>','\rm');
TString1 = strrep(TString1,'_','\_');
TString2 = strrep(TString2,'<strong>','\bf');
TString2 = strrep(TString2,'</strong>','\rm');
TString2 = strrep(TString2,'_','\_');
TString3 = strrep(TString3,'<strong>','\bf');
TString3 = strrep(TString3,'</strong>','\rm');
TString3 = strrep(TString3,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
annotation(gcf,'Textbox','String',TString1,'Interpreter','Tex',...
    'FontName',FixedWidth,'LineStyle','none','Fontsize',15,'Units','Normalized','Position',[0.62 -0.1 0.7 0.5]);
annotation(gcf,'Textbox','String',TString2,'Interpreter','Tex',...
    'FontName',FixedWidth,'LineStyle','none','Fontsize',15,'Units','Normalized','Position',[0.62 -0.22 0.7 0.5]);
annotation(gcf,'Textbox','String',TString3,'Interpreter','Tex',...
    'FontName',FixedWidth,'LineStyle','none','Fontsize',15,'Units','Normalized','Position',[0.62 -0.34 0.7 0.5]);

%% Stability region for the Radau5 method
gamma = (2-sqrt(2))/2;
A = [[0,0,0];...
     [gamma,gamma,0];...
     [sqrt(2)/4,sqrt(2)/4,gamma]];
b = [sqrt(2)/4;sqrt(2)/4;gamma];
e = [1;1;1];

alpha = -20:0.1:20;
beta = -20:0.1:20;
nreal = length(alpha);
nimag = length(beta);

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
colorbar;
% contourf(alpha,beta,absR,[0,1])
% colormap(map)
plot([0,0],[-20,20],'k-')
plot([-20,20],[0,0],'k-')
axis image;
axis xy;
xlabel('real');
ylabel('imag');
title('The Stability region (|R(z)|) for the ESDIRK23')
set(gca,'fontsize',fs)

