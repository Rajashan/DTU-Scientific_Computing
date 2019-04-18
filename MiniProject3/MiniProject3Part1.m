%% MiniProject3 Part 1
clear all; close all; clc

%% Test Equation

h = 0.1;
tspan = linspace(0,5,10^4);
x0 = 1;
lambda = -1;
s = 1;
fun = @(t,x) lambda(s)*x;

[Tout,Xout,E] = RungeKutta32Solver(fun,tspan,x0,h);


%%
clear variables;

h = 0.1;
tspan = linspace(0,5,10^4);
x0 = 1;
lambda = -1;
s = 1;
fun = @(t,x) lambda(s)*x;

x_ode45 = ode45(fun,tspan,x0);
x_ode15s = ode15s(fun,tspan,x0);

[Tout,Xout] = RungeKutta32Solver(fun,tspan,x0,h);
[Xa] = x0*exp(lambda(s)*(Tout-Tout(1)));

figure(1)
subplot(2,1,1)
plot(Tout,Xout,'Linewidth',3,'Linestyle','-','Color','Blue')
hold on
grid on
plot(x_ode45.x,x_ode45.y,'Linewidth',3,'Linestyle','-.','Color','Red')
plot(x_ode15s.x,x_ode15s.y,'Linewidth',3,'Linestyle','--','Color','Black')
plot(Tout,Xa,'Linewidth',3,'Linestyle','--','Color','Green')
title({['\textbf{Test Equation:} $\lambda$ = ',num2str(lambda(s))],[]},'Interpreter','Latex','Fontsize',18)
lg = legend(['ERK3(2) Solution (Steps = ',num2str(length(Tout)),')'],['ODE45 Solution (Steps = ',num2str(length(x_ode45.x)),')'],['ODE15s Solution (Steps = ',num2str(length(x_ode15s.x)),')'],'Analytical Solution','location','northeast');
ylim([0 max(Xout)])
xlim([0 5])
set(lg,'Interpreter','Latex','FontSize',20)
xlabel('$t$','Interpreter','Latex')
ylabel('$x(t)$','Interpreter','Latex')
set(gca, 'FontSize', 20)

subplot(2,1,2)
plot(Tout,Xout,'Linewidth',3,'Linestyle','-','Color','Blue')
hold on
grid on
plot(x_ode45.x,x_ode45.y,'Linewidth',3,'Linestyle','-.','Color','Red')
plot(x_ode15s.x,x_ode15s.y,'Linewidth',3,'Linestyle','--','Color','Black')
plot(Tout,Xa,'Linewidth',3,'Linestyle','--','Color','Green')
title({['\textbf{Test Equation:} $\lambda$ = ',num2str(lambda(s))],[]},'Interpreter','Latex','Fontsize',18)
lg = legend(['ERK3(2) Solution (Steps = ',num2str(length(Tout)),')'],['ODE45 Solution (Steps = ',num2str(length(x_ode45.x)),')'],['ODE15s Solution (Steps = ',num2str(length(x_ode15s.x)),')'],'Analytical Solution','location','northeast');
ylim([0 max(Xout)])
xlim([0 5])
set(lg,'Interpreter','Latex','FontSize',20)
xlabel('$t$','Interpreter','Latex')
ylabel('$x(t)$','Interpreter','Latex')
set(gca, 'FontSize', 20)
ylim([0.48 0.51])
xlim([0.7225 0.7226])

set(gcf,'units','points','position',[150,0,1500,1000])
print('TestEqRK32','-depsc')

%% Convergence Rate
clear variables;

h = (1/2).^(1:7+2);
tspan = linspace(0,5,10^4);
x0 = 1;
lambda = -1;
fun = @(t,x) lambda*x;

for s = 1 : length(h)

[Tout,Xout,E] = RungeKutta32Solver(fun,tspan,x0,h(s));
%Calculate LTE after second point
LTE(s) = E(2);
end
 
figure(1)
%Local Error
Plot1 = plot(h,LTE,'linewidth',4,'Color','Blue');
hold on;grid on;
Plot2 = scatter(h,LTE,100,'filled','MarkerFaceColor','Blue');
set(get(get(Plot2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
Plot3 = plot(h,h.^3,'linewidth',4,'Color','Blue','Linestyle','--');
set(gca,'XScale','log','YScale','log')
title({'\textbf{Convergence Plot for Explicit Runge-Kutta3(2) Method}',['\textbf{Test Equation:}','$\lambda$ = ',num2str(lambda)]},'Interpreter','Latex')
xlabel('Step Size $h$','Interpreter','Latex')
ylabel('Error $|e(x(t))|$','Interpreter','Latex')
lg = legend('Local Truncation Error','Reference Line $f(h) = h^{3}$');
set(lg,'Interpreter','Latex','FontSize',20,'location','southeast')
set(gca, 'FontSize', 17)
set(gcf,'units','points','position',[450,150,800,500])
print('TestEqERK32_Convergence','-depsc')

%% Stability for Test Equation
clear variables;
%Butcher Tableau Information for ERK4C
s = 3;
A = [0 0 0 ; 1/2 0 0 ; -1 2 0];
b = [1/6 ; 2/3 ; 1/6];
c = [ 0 ; 1/2 ; 1];

%Interval of Complex Region
alpha = linspace(-4,4,1000);
beta = linspace(-4,4,1000);
nreal = length(alpha);
nimag = length(beta);

%Initial Values
I = eye(size(A));
e = ones(size(A,1),1);
absR = zeros(length(nimag),length(nreal));
%Compute the transfer functions R(z)
for kreal = 1:nreal
   for kimag = 1:nimag
        z = alpha(kreal) + 1*i*beta(kimag);
        tmp = (I-z*A)\e;
        R = 1 + z*b'*tmp;
        absR(kimag,kreal) = abs(R);
   end
end
%%
%Plot Transfer Function
figure(2)
fs = 14;
imagesc(alpha,beta,absR,[0 1]);
grid on
colorbar
axis image
axis xy
xlabel('$Re(z)$','Interpreter','Latex','fontsize',fs);
ylabel('$Im(z)$','Interpreter','Latex','fontsize',fs);
title('$|R(z)|$','Interpreter','Latex','fontsize',fs)
set(gca, 'FontSize', 18)
set(gcf,'units','points','position',[450,150,600,500])
print('TestEqStab_ERK32','-depsc')

%% Van Der Pol Problem
clear variables;
mu = 3;
x0 = [0.5 0.5]';
tspan = linspace(0,50);
h = 0.01;
options = odeset('RelTol',10^(-6),'AbsTol',10^(-6));

for s = 1:length(mu)
fun = @(t,x) [x(2) ; mu(s) * ( 1 - x(1)^2) * x(2) - x(1)];
figure(s)
[Tout,Xout] = RungeKutta32Solver(@VanDerPolJac,tspan,x0,h,mu(s));
%x_ode45 = ode45(fun,tspan,x0,options);
x_ode15s = ode45(fun,tspan,x0,options);

subplot(2,2,1)
plot(Tout,Xout(:,1),'Linestyle','-','Color','Blue','LineWidth',3)
hold on
grid on
%plot(x_ode45.x,x_ode45.y(1,:),'Linestyle','--','Color','Red','LineWidth',3)
plot(x_ode15s.x,x_ode15s.y(1,:),'Linestyle','--','Color','Green','LineWidth',1)
title({'\textbf{Van Der Pol} ',['$x_{1}(t) \quad \mu$ = ',num2str(mu(s))]},'Interpreter','Latex','FontSize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$x_{1}(t)$','Interpreter','Latex')
set(gca, 'FontSize', 18)

subplot(2,2,2)
plot(Tout,Xout(:,2),'Linestyle','-','Color','Blue','LineWidth',3)
hold on
grid on
%plot(x_ode45.x,x_ode45.y(2,:),'Linestyle','--','Color','Red','LineWidth',3)
plot(x_ode15s.x,x_ode15s.y(2,:),'Linestyle','--','Color','Green','LineWidth',1)
title({'\textbf{Van Der Pol} ',['$x_{2}(t) \quad \mu$ = ',num2str(mu(s))]},'Interpreter','Latex','FontSize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$x_{2}(t)$','Interpreter','Latex')
set(gca, 'FontSize', 18)

subplot(2,2,3:4)
plot(Xout(:,1),Xout(:,2),'Linestyle','-','Color','Blue','LineWidth',3)
hold on
grid on
%plot(x_ode45.y(1,:),x_ode45.y(2,:),'Linestyle','-','Color','Red','LineWidth',3)
plot(x_ode15s.y(1,:),x_ode15s.y(2,:),'Linestyle','-','Color','Green','LineWidth',1)
title({'\textbf{Van Der Pol Phase Plot} ',['$\mu$ = ',num2str(mu(s))]},'Interpreter','Latex','FontSize',18)
xlabel('$x_{1}(t)$','Interpreter','Latex')
ylabel('$x_{2}(t)$','Interpreter','Latex')
lg = legend(['ERK3(2) ','(Steps = ',num2str(length(Xout)),')'],['ODE15s ','(Steps = ',num2str(length(x_ode15s.y)),')']);
set(lg,'Interpreter','Latex','FontSize',20,'location','southeast')
set(gca, 'FontSize', 18)


set(gcf,'units','points','position',[450,150,1000,700])
print(['VanDerPol_RK32',num2str(s)],'-depsc')

end

%% Van Der Pol Problem
clear all;

reltol = [5*10^-4];
abstol = [5*10^-4];
x0 = [0.5 0.5]';
tspan = linspace(0,50);
h = 0.1;

mu = [3 10];
%options = [];
options = odeset('RelTol',10^-6,'AbsTol',10^-6);

for s = 1:1
fun = @(t,x) [x(2) ; mu(s) * ( 1 - x(1)^2) * x(2) - x(1)];
figure(s)
[Tout,Xout,H,tau,R] = RungeKutta32AdaptiveSolver(@VanDerPolJac,tspan,x0,h,abstol(s),reltol(s),mu(s));
%x_ode45 = ode45(fun,tspan,x0,options);
x_ode15s = ode45(fun,tspan,x0,options);

subplot(2,2,1)
plot(Tout,Xout(1,:),'Linestyle','-','Color','Blue','LineWidth',3)
hold on
grid on
%plot(x_ode45.x,x_ode45.y(1,:),'Linestyle','--','Color','Red','LineWidth',3)
plot(x_ode15s.x,x_ode15s.y(1,:),'Linestyle','--','Color','Green','LineWidth',1)
title({'\textbf{Van Der Pol} ',['$x_{1}(t) \quad \mu$ = ',num2str(mu(s))]},'Interpreter','Latex','FontSize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$x_{1}(t)$','Interpreter','Latex')
set(gca, 'FontSize', 18)

subplot(2,2,2)
plot(Tout,Xout(2,:),'Linestyle','-','Color','Blue','LineWidth',3)
hold on
grid on
%plot(x_ode45.x,x_ode45.y(2,:),'Linestyle','--','Color','Red','LineWidth',3)
plot(x_ode15s.x,x_ode15s.y(2,:),'Linestyle','--','Color','Green','LineWidth',1)
title({'\textbf{Van Der Pol} ',['$x_{2}(t) \quad \mu$ = ',num2str(mu(s))]},'Interpreter','Latex','FontSize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$x_{2}(t)$','Interpreter','Latex')
set(gca, 'FontSize', 18)

subplot(2,2,3:4)
plot(Xout(1,:),Xout(2,:),'Linestyle','-','Color','Blue','LineWidth',3)
hold on
grid on
%plot(x_ode45.y(1,:),x_ode45.y(2,:),'Linestyle','-','Color','Red','LineWidth',3)
plot(x_ode15s.y(1,:),x_ode15s.y(2,:),'Linestyle','-','Color','Green','LineWidth',1)
title({'\textbf{Van Der Pol Phase Plot} ',['$\mu$ = ',num2str(mu(s))]},'Interpreter','Latex','FontSize',18)
xlabel('$x_{1}(t)$','Interpreter','Latex')
ylabel('$x_{2}(t)$','Interpreter','Latex')
lg = legend(['ERK3(2) ','(Steps = ',num2str(length(Xout)),')'],['ODE15s ','(Steps = ',num2str(length(x_ode15s.y)),')']);
set(lg,'Interpreter','Latex','FontSize',20,'location','southeast')
set(gca, 'FontSize', 18)


set(gcf,'units','points','position',[450,150,1000,700])
print(['VanDerPolERK32',num2str(s)],'-depsc')
end

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
print('VanDerPolERK32_HR','-depsc')