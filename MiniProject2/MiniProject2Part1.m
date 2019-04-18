%% MiniProject 2 Part 1
clear variables; close all; clc;

%% Test Equation with Classical Runge-Kutta
clear variables;

h = 0.1;
tspan = linspace(0,3,10^4);
x0 = 1;
lambda = [-1 -10 -20 -25 1 10 20 25];
%lambda = [-1 -1];

for s = 1:length(lambda)
fun = @(t,x) lambda(s)*x;
x_ode45 = ode45(fun,tspan,x0);
x_ode15s = ode15s(fun,tspan,x0);
[Tout,Xout] = ClassicalRungeKuttaSolver(fun,tspan,x0,h);
[Xa] = x0*exp(lambda(s)*(Tout-Tout(1)));

subplot(2,length(lambda)/2,s)
%figure(1)
plot(Tout,Xout,'Linewidth',3,'Linestyle','-','Color','Blue')
hold on
grid on
plot(x_ode45.x,x_ode45.y,'Linewidth',3,'Linestyle','-.','Color','Red')
%plot(x_ode15s.x,x_ode15s.y,'Linewidth',3,'Linestyle','--','Color','Black')
plot(Tout,Xa,'Linewidth',3,'Linestyle','--','Color','Green')
title({['\textbf{Test Equation:} $\lambda$ = ',num2str(lambda(s))],[]},'Interpreter','Latex','Fontsize',18)
if s > 4.5
    lg = legend('ERK4C Solution','ODE45 Solution','Analytical Solution','location','northwest');
    ylim([0 max(Xout)])
    xlim([0 5])
else
    lg = legend('ERK4C Solution','ODE45 Solution','Analytical Solution','location','northeast');
    ylim([0 max(Xout)])
    xlim([0 5])
end
set(lg,'Interpreter','Latex','FontSize',14)
xlabel('$t$','Interpreter','Latex')
ylabel('$x(t)$','Interpreter','Latex')
set(gca, 'FontSize', 15)

end
set(gcf,'units','points','position',[150,0,1500,1000])
print('TestEq','-depsc')

%% Convergence Rate
clear variables;

h = (1/2).^(1:7+2);
tspan = linspace(0,2);
x0 = 1;
p = [1 -1];
lambda = p(1);
fun = @(t,x) lambda*x;

for s = 1 : length(h)

[Tout,Xout] = ClassicalRungeKuttaSolver(fun,tspan,x0,h(s));
Xa = x0*exp(lambda*(Tout-Tout(1)));

%Calculate LTE at first point
LTE(s) = abs(Xa(2) - Xout(2));
%LTE(s) = abs( Xout(round(end/2)) - Xa(round(end/2)) );

%Calculate GTE
GTE(s) = abs(Xout(end) - Xa(end));
%GTE(s) = sum(abs(Xa - Xout));
end
 
figure(1)
%Local Error
Plot1 = plot(h,LTE,'linewidth',4,'Color','Blue');
hold on;grid on;
Plot2 = scatter(h,LTE,100,'filled','MarkerFaceColor','Blue');
set(get(get(Plot2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
Plot3 = plot(h,h.^5/20,'linewidth',4,'Color','Blue','Linestyle','--');
%Global Error
Plot4 = plot(h,GTE*10^2,'linewidth',4,'Color','Red');
Plot5 = scatter(h,GTE*10^2,100,'filled','MarkerFaceColor','Red');
set(get(get(Plot5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
Plot6 = plot(h,h.^4,'linewidth',4,'Color','Red','Linestyle','--');
%Options, Legends, Titles etc.
set(gca,'XScale','log','YScale','log')
title({'\textbf{Convergence Plot for Classical Runge-Kutta Method}',['$\lambda$ = ',num2str(lambda)]},'Interpreter','Latex')
xlabel('Step Size $h$','Interpreter','Latex')
ylabel('Error $|e(x(t))|$','Interpreter','Latex')
lg = legend('Local Truncation Error','Reference Line $f(h) = h^{5}$','Global Truncation Error','Reference Line $f(h) = h^{4}$');
set(lg,'Interpreter','Latex','FontSize',20,'location','southeast')
xlim([10^-3 10^0])
ylim([10^(-18) 10^1]);
set(gca, 'FontSize', 17)
set(gcf,'units','points','position',[450,150,800,500])
print('TestEqConvergence','-depsc')

%%
fun = @(t,x) lambda*x;
lambda = p(2);
for s = 1 : length(h)
tic
[Tout,Xout] = ClassicalRungeKuttaSolver(fun,tspan,x0,h(s));
Xa = x0*exp(lambda*(Tout-Tout(1)));

%Calculate LTE at first point
LTE(s) = abs(Xout(2) - Xa(2));
%Calculate GTE at end point
GTE(s) = abs(Xout(end) - Xa(end));
toc
end

figure(2)
%Local Error
Plot1 = plot(h,LTE/10^5,'linewidth',3,'Color','Blue');
hold on;grid on;
Plot2 = scatter(h,LTE/10^5,'filled','MarkerFaceColor','Blue');
set(get(get(Plot2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
Plot3 = plot(h,h.^5/10^8,'linewidth',3);
%Global Error
Plot4 = plot(h,GTE*10,'linewidth',3,'Color','Red');
Plot5 = scatter(h,GTE*10,'filled','MarkerFaceColor','Red');
set(get(get(Plot5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
Plot6 = plot(h,h.^4/10^2,'linewidth',3);
%Options, Legends, Titles etc.
set(gca,'XScale','log','YScale','log')
title({'Convergence Plot for Classical Runge-Kutta Method',['$\lambda$ = ',num2str(lambda)]},'Interpreter','Latex','FontSize',18)
xlabel('Step Size h')
ylabel('Error $e(x_{n})$','Interpreter','Latex')
lg = legend('Local Truncation Error','Reference Line $h = h^{5}$','Global Truncation Error','Reference Line $h = h^{4}$');
set(lg,'Interpreter','Latex','FontSize',15,'location','southeast')
%xlim([10^-3 10^1])
%ylim([10^(-16) 10^0]);
set(gca, 'FontSize', 15)

%% Stability for Test Equation
clear variables;
%Butcher Tableau Information for ERK4C
s  = 4;
A = [0 0 0 0 ; 1/2 0 0 0 ; 0 1/2 0 0 ; 0 0 1 0];
b  = [1/6 ; 1/3 ; 1/3 ; 1/6];
c  = [0; 1/2 ; 1/2 ; 1];

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
print('TestEqStab','-depsc')

%% Van Der Pol Problem
clear variables;
mu = [3 10];
%mu = 3;
x0 = [0.5 0.5]';
tspan = linspace(0,50);
h = 0.01;
options = odeset('RelTol',10^(-6),'AbsTol',10^(-6));
%options = [];
for s = 1:length(mu)
fun = @(t,x) [x(2) ; mu(s) * ( 1 - x(1)^2) * x(2) - x(1)];
figure(s)
[Tout,Xout] = ClassicalRungeKuttaSolver(@VanDerPolJac,tspan,x0,h,mu(s));
x_ode45 = ode45(fun,tspan,x0,options);
x_ode15s = ode45(fun,tspan,x0,options);

subplot(2,2,1)
plot(Tout,Xout(:,1),'Linestyle','-','Color','Blue','LineWidth',3)
hold on
grid on
plot(x_ode45.x,x_ode45.y(1,:),'Linestyle','--','Color','Red','LineWidth',3)
plot(x_ode15s.x,x_ode15s.y(1,:),'Linestyle','--','Color','Green','LineWidth',1)
title({'\textbf{Van Der Pol} ',['$x_{1}(t) \quad \mu$ = ',num2str(mu(s))]},'Interpreter','Latex','FontSize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$x_{1}(t)$','Interpreter','Latex')
set(gca, 'FontSize', 18)

subplot(2,2,2)
plot(Tout,Xout(:,2),'Linestyle','-','Color','Blue','LineWidth',3)
hold on
grid on
plot(x_ode45.x,x_ode45.y(2,:),'Linestyle','--','Color','Red','LineWidth',3)
plot(x_ode15s.x,x_ode15s.y(2,:),'Linestyle','--','Color','Green','LineWidth',1)
title({'\textbf{Van Der Pol} ',['$x_{2}(t) \quad \mu$ = ',num2str(mu(s))]},'Interpreter','Latex','FontSize',18)
xlabel('$t$','Interpreter','Latex')
ylabel('$x_{2}(t)$','Interpreter','Latex')
set(gca, 'FontSize', 18)

subplot(2,2,3:4)
plot(Xout(:,1),Xout(:,2),'Linestyle','-','Color','Blue','LineWidth',3)
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

%% Prey Predator
clear variables;
x0 = [3 1];
p = [2 2];
tspan = linspace(0,15);
h = 0.01;
options = odeset('RelTol',10^(-9),'AbsTol',10^(-9));
%options = [];
for s = 1:1
fun = @(t,x) [p(1) * (1 - x(2)) * x(1) ; -p(2) * (1 - x(1)) * x(2)];
figure(s)
[Tout,Xout] = ClassicalRungeKuttaSolver(@preyPredator,tspan,x0,h,p);
x_ode45 = ode45(fun,tspan,x0,options);
x_ode15s = ode45(fun,tspan,x0,options);

subplot(2,2,1)
plot(Tout,Xout(:,1),'Linestyle','-','Color','Blue','LineWidth',3)
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
plot(Tout,Xout(:,2),'Linestyle','-','Color','Blue','LineWidth',3)
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
plot(Xout(:,1),Xout(:,2),'Linestyle','-','Color','Blue','LineWidth',3)
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
xlim([min(Xout(:,1)),max(Xout(:,1))])


set(gcf,'units','points','position',[450,150,1000,700])
print('PreyPredator','-depsc')

end