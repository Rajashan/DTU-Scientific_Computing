%% MiniProject 2 Part 2
clear variables; close all; clc;

%% Test Equation with Adaptive Runge-Kutta
clear variables;
%reltol = 10^(-3);
%abstol = 10^(-3);
%reltol = 0.75;
%abstol = 0.75;
reltol = 0.1;
abstol = 0.1;

h0 = 0.1;
tspan = linspace(0,5,10^4);
x0 = 1;
%lambda = [1 10 15 25];
lambda = [25 25];
options = odeset('RelTol',10^-3,'AbsTol',10^-3);
%options = [];
for s = 1:length(lambda)
fun = @(t,x) lambda(s)*x;
x_ode45 = ode45(fun,tspan,x0,options);
%x_ode15s = ode15s(fun,tspan,x0);
[Tout,Xout,H,tau,R] = AdaptiveRungeKuttaSolver(fun, tspan, x0, h0, abstol, reltol);
[Xa] = x0*exp(lambda(s)*(tspan-tspan(1)));

subplot(2,length(lambda)/2,s)
%figure(1)
plot(Tout,Xout,'Linewidth',3,'Linestyle','-','Color','Blue')
hold on
grid on
plot(x_ode45.x,x_ode45.y,'Linewidth',3,'Linestyle','-.','Color','Red')
%plot(x_ode15s.x,x_ode15s.y,'Linewidth',3,'Linestyle','--','Color','Black')
plot(tspan,Xa,'Linewidth',3,'Linestyle','--','Color','Green')
title({['\textbf{Test Equation:} $\lambda$ = ',num2str(lambda(s))],[]},'Interpreter','Latex','Fontsize',18)
%if s > 4.5
    lg = legend(['ERK4C Solution',' (Steps = ',num2str(length(Tout)),')'],['ODE45 Solution',' (Steps = ',num2str(length(x_ode45.x)),')'],'Analytical Solution','location','northwest');
    ylim([0 max(Xout)])
    %xlim([0 0.2])
%else
 %   lg = legend('ERK4C Solution','Analytical Solution','ODE45 Solution','location','northeast');
 %  ylim([0 max(Xout)])
 %  xlim([0 5.5])
%end
set(lg,'Interpreter','Latex','FontSize',14)
xlabel('$t$','Interpreter','Latex')
ylabel('$x(t)$','Interpreter','Latex')
set(gca, 'FontSize', 15)
s
length(Tout)
end
set(gcf,'units','points','position',[150,0,1500,1000])
print('TestEqAdaptive','-depsc')
%%
figure(2)
subplot(2,2,1)
plot(Tout,Xout,'Linewidth',5,'Linestyle','-','Color','Blue')
%scatter(Tout,abs(tau),10,'MarkerFaceColor','Blue')
hold on
grid on
title(['\textbf{Test Equation:} $\lambda$ = ',num2str(lambda(end))],'Interpreter','Latex','Fontsize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$x(t_{k})$','Interpreter','Latex')
set(gca, 'FontSize', 15)

subplot(2,2,3)
%plot(Tout,Xout,'Linewidth',3,'Linestyle','-','Color','Blue')
scatter(Tout,abs(tau),50,'MarkerFaceColor','Blue')
hold on
grid on
title(['\textbf{Error E:} $\lambda$ = ',num2str(lambda(end))],'Interpreter','Latex','Fontsize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$|E(t_{k})|$','Interpreter','Latex')
set(gca, 'FontSize', 15)

subplot(2,2,2)
%plot(Tout,H,'Linewidth',3,'Linestyle','-','Color','Red') 
scatter(Tout,H,50,'MarkerFaceColor','Blue')
hold on
grid on
title(['\textbf{Step Size H:} $\lambda$ = ',num2str(lambda(end))],'Interpreter','Latex','Fontsize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$h_{k}$','Interpreter','Latex')
set(gca, 'FontSize', 15)

subplot(2,2,4)
%plot(Tout,H,'Linewidth',3,'Linestyle','-','Color','Red') 
scatter(Tout,R,50,'MarkerFaceColor','Blue')
hold on
grid on
title(['\textbf{Step Acceptor R:} $\lambda$ = ',num2str(lambda(end))],'Interpreter','Latex','Fontsize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$R_{k}$','Interpreter','Latex')
set(gca, 'FontSize', 15)
set(gcf,'units','points','position',[150,0,1500,1000])
print('TestEqAdaptiveHR','-depsc')

%% Van Der Pol Problem
clear all;

reltol = [7*10^(-2), 6*10^(-2)];
abstol = [7*10^(-2), 6*10^(-2)];
x0 = [0.5 0.5]';
tspan = linspace(0,50);
h = 0.1;

mu = [3 10];
options = [];

for s = 1:length(mu)
fun = @(t,x) [x(2) ; mu(s) * ( 1 - x(1)^2) * x(2) - x(1)];
figure(s)
[Tout,Xout,H,tau,R] = AdaptiveRungeKuttaSolver(@VanDerPolJac,tspan,x0,h,abstol(s),reltol(s),mu(s));
%options = odeset('RelTol',10^-6,'AbsTol',10^-6);
options = odeset('RelTol',10^-6,'AbsTol',10^-6);
x_ode45 = ode45(fun,tspan,x0,options);
x_ode15s = ode45(fun,tspan,x0,options);

subplot(2,2,1)
plot(Tout,Xout(1,:),'Linestyle','-','Color','Blue','LineWidth',3)
hold on
grid on
plot(x_ode45.x,x_ode45.y(1,:),'Linestyle','--','Color','Red','LineWidth',3)
plot(x_ode15s.x,x_ode15s.y(1,:),'Linestyle','--','Color','Green','LineWidth',1)
title({'\textbf{Van Der Pol} ',['$x_{1}(t) \quad \mu$ = ',num2str(mu(s))]},'Interpreter','Latex','FontSize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$x_{1}(t)$','Interpreter','Latex')
set(gca, 'FontSize', 18)

subplot(2,2,2)
plot(Tout,Xout(2,:),'Linestyle','-','Color','Blue','LineWidth',3)
hold on
grid on
plot(x_ode45.x,x_ode45.y(2,:),'Linestyle','--','Color','Red','LineWidth',3)
plot(x_ode15s.x,x_ode15s.y(2,:),'Linestyle','--','Color','Green','LineWidth',1)
title({'\textbf{Van Der Pol} ',['$x_{2}(t) \quad \mu$ = ',num2str(mu(s))]},'Interpreter','Latex','FontSize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$x_{2}(t)$','Interpreter','Latex')
set(gca, 'FontSize', 18)

subplot(2,2,3:4)
plot(Xout(1,:),Xout(2,:),'Linestyle','-','Color','Blue','LineWidth',3)
hold on
grid on
plot(x_ode45.y(1,:),x_ode45.y(2,:),'Linestyle','-','Color','Red','LineWidth',3)
plot(x_ode15s.y(1,:),x_ode15s.y(2,:),'Linestyle','-','Color','Green','LineWidth',1)
title({'\textbf{Van Der Pol Phase Plot} ',['$\mu$ = ',num2str(mu(s))]},'Interpreter','Latex','FontSize',18)
xlabel('$x_{1}(t)$','Interpreter','Latex')
ylabel('$x_{2}(t)$','Interpreter','Latex')
lg = legend(['ERK4C ','(Steps = ',num2str(length(Xout)),')'],['ODE45 ','(Steps = ',num2str(length(x_ode45.y)),')'],['ODE15s ','(Steps = ',num2str(length(x_ode15s.y)),')']);
set(lg,'Interpreter','Latex','FontSize',20,'location','southeast')
set(gca, 'FontSize', 18)


set(gcf,'units','points','position',[450,150,1000,700])
print(['VanDerPol',num2str(s)],'-depsc')

end
%%
figure(3)
subplot(2,2,1)
plot(Tout,Xout(2,:),'Linewidth',5,'Linestyle','-','Color','Blue')
%scatter(Tout,abs(tau),10,'MarkerFaceColor','Blue')
hold on
grid on
title({['\textbf{Van Der Pol:} $\mu$ = ',num2str(mu(end))],'$x_{1}(t)$'},'Interpreter','Latex','Fontsize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$x_{1}(t_{k})$','Interpreter','Latex')
set(gca, 'FontSize', 18)

subplot(2,2,3)
plot(Tout,abs(tau),'Linewidth',3,'Linestyle','-','Color','Blue')
%scatter(Tout,abs(tau),50,'MarkerFaceColor','Blue')
hold on
grid on
title(['\textbf{Error E:} $\mu$ = ',num2str(mu(end))],'Interpreter','Latex','Fontsize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$|E(t_{k})|$','Interpreter','Latex')
set(gca, 'FontSize', 18)

subplot(2,2,2)
plot(Tout,H,'Linewidth',3,'Linestyle','-','Color','Blue') 
%scatter(Tout,H,50,'MarkerFaceColor','Blue')
hold on
grid on
title(['\textbf{Step Size H:} $\mu$ = ',num2str(mu(end))],'Interpreter','Latex','Fontsize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$h_{k}$','Interpreter','Latex')
set(gca, 'FontSize', 18)

subplot(2,2,4)
plot(Tout,R,'Linewidth',2,'Linestyle','-','Color','Blue') 
%scatter(Tout,R,50,'MarkerFaceColor','Blue')
hold on
grid on
title(['\textbf{Step Acceptor R:} $\mu$ = ',num2str(mu(end))],'Interpreter','Latex','Fontsize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$R_{k}$','Interpreter','Latex')
set(gca, 'FontSize', 18)
set(gcf,'units','points','position',[150,0,1500,1000])
print('VanDerPolAdaptiveHR','-depsc')

%% Prey Predator
clear variables; close all;
reltol = 3*10^(-2);
abstol = 3*10^(-2);

x0 = [3 1];
p = [2 2];
tspan = linspace(0,15);
h = 0.1;
options = odeset('RelTol',10^(-9),'AbsTol',10^(-9));
%options = [];
for s = 1:1
fun = @(t,x) [p(1) * (1 - x(2)) * x(1) ; -p(2) * (1 - x(1)) * x(2)];
figure(s)
[Tout,Xout,H,tau,R] = AdaptiveRungeKuttaSolver(@preyPredator,tspan,x0,h,abstol(s),reltol(s),p);
x_ode45 = ode45(fun,tspan,x0,options);
x_ode15s = ode45(fun,tspan,x0,options);

subplot(2,2,1)
plot(Tout,Xout(1,:),'Linestyle','-','Color','Blue','LineWidth',3)
hold on
grid on
plot(x_ode45.x,x_ode45.y(1,:),'Linestyle','--','Color','Red','LineWidth',3)
plot(x_ode15s.x,x_ode15s.y(1,:),'Linestyle','--','Color','Green','LineWidth',1)
title({'\textbf{Prey Predator} ',['$x_{1}(t) \quad p$ =  [',num2str(p),']']},'Interpreter','Latex','FontSize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$x_{1}(t)$','Interpreter','Latex')
set(gca, 'FontSize', 18)
xlim([tspan(1), tspan(end)])

subplot(2,2,2)
plot(Tout,Xout(2,:),'Linestyle','-','Color','Blue','LineWidth',3)
hold on
grid on
plot(x_ode45.x,x_ode45.y(2,:),'Linestyle','--','Color','Red','LineWidth',3)
plot(x_ode15s.x,x_ode15s.y(2,:),'Linestyle','--','Color','Green','LineWidth',1)
title({'\textbf{Prey Predator} ',['$x_{2}(t) \quad p$ =  [',num2str(p),']']},'Interpreter','Latex','FontSize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$x_{2}(t)$','Interpreter','Latex')
set(gca, 'FontSize', 18)
xlim([tspan(1), tspan(end)])

subplot(2,2,3:4)
plot(Xout(1,:),Xout(2,:),'Linestyle','-','Color','Blue','LineWidth',3)
hold on
grid on
plot(x_ode45.y(1,:),x_ode45.y(2,:),'Linestyle','-','Color','Red','LineWidth',3)
plot(x_ode15s.y(1,:),x_ode15s.y(2,:),'Linestyle','-','Color','Green','LineWidth',1)
title({'\textbf{Prey Predator Phase Plot} ',['$p$ = [',num2str(p),']']},'Interpreter','Latex','FontSize',18)
xlabel('$x_{1}(t)$','Interpreter','Latex')
ylabel('$x_{2}(t)$','Interpreter','Latex')
lg = legend(['ERK4C ','(Steps = ',num2str(length(Xout)),')'],['ODE45 ','(Steps = ',num2str(length(x_ode45.y)),')'],['ODE15s ','(Steps = ',num2str(length(x_ode15s.y)),')']);
set(lg,'Interpreter','Latex','FontSize',20,'location','northeast')
set(gca, 'FontSize', 18)
xlim([min(Xout(1,:)),max(Xout(2,:))])

set(gcf,'units','points','position',[450,150,1000,700])
print('PreyPredator','-depsc')

end

%%
figure()
subplot(2,2,1)
plot(Tout,Xout(2,:),'Linewidth',5,'Linestyle','-','Color','Blue')
%scatter(Tout,abs(tau),10,'MarkerFaceColor','Blue')
hold on
grid on
title({['\textbf{Prey Predator:} ','$p$ = [',num2str(p),']'],'\qquad $x_{1}(t)$'},'Interpreter','Latex','FontSize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$x_{1}(t_{k})$','Interpreter','Latex')
set(gca, 'FontSize', 18)

subplot(2,2,3)
plot(Tout,abs(tau),'Linewidth',3,'Linestyle','-','Color','Blue')
%scatter(Tout,abs(tau),50,'MarkerFaceColor','Blue')
hold on
grid on
title(['\textbf{Error E:} ','$p$ = [',num2str(p),']'],'Interpreter','Latex','Fontsize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$|E(t_{k})|$','Interpreter','Latex')
set(gca, 'FontSize', 18)

subplot(2,2,2)
plot(Tout,H,'Linewidth',3,'Linestyle','-','Color','Blue') 
%scatter(Tout,H,50,'MarkerFaceColor','Blue')
hold on
grid on
title(['\textbf{Step Size H:} ','$p$ = [',num2str(p),']'],'Interpreter','Latex','Fontsize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$h_{k}$','Interpreter','Latex')
set(gca, 'FontSize', 18)

subplot(2,2,4)
plot(Tout,R,'Linewidth',3,'Linestyle','-','Color','Blue') 
%scatter(Tout,R,50,'MarkerFaceColor','Blue')
hold on
grid on
title(['\textbf{Step Accepter R:} ','$p$ = [',num2str(p),']'],'Interpreter','Latex','Fontsize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$R_{k}$','Interpreter','Latex')
set(gca, 'FontSize', 18)
set(gcf,'units','points','position',[150,0,1500,1000])
print('PreyPredatorAdaptiveHR','-depsc')
