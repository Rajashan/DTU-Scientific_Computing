%% Solve Van Der Pol using Explicit Euler Fixed Time Step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables;close all;clc

ta = 0;
tb = 25;

%Numerical Solution
N = 10^3;
x0 = [0.5 0.5];
mu = 3;
tspan = [ta tb];
[T,X] = ExplicitEuler(@VanDerPolJac,tspan,N,x0,mu);

%Plotting
figure(1)
subplot(2,2,1)
plot(T,X(1,:),'b-','linewidth',2)
xlabel('t')
ylabel('$x_{1}$','Interpreter','Latex','Rotation',0,'FontSize',15)
title('$x_{1}(t)$','Interpreter','Latex','FontSize',18)

subplot(2,2,3)
plot(T,X(2,:),'b-','linewidth',2)
xlabel('t')
ylabel('$x_{2}$','Interpreter','Latex','Rotation',0,'FontSize',15)
title('$x_{2}(t)$','Interpreter','Latex','FontSize',18)

subplot(2,2,[2 4])
plot(X(1,:),X(2,:),'b-','linewidth',2)
xlabel('$x_{1}$','Interpreter','Latex','Rotation',0,'FontSize',15)
ylabel('$x_{2}$','Interpreter','Latex','Rotation',0,'FontSize',15)
title('Phase Plot','Interpreter','Latex','FontSize',18)

%% Solve Van Der Pol using Explicit Euler Adaptive Time Step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Numerical Solution
tspan = [0 25];
x0 = [0.5 0.5];
h0 = 10^-2;
abstol = 10^(-5);
reltol = 10^(-5);
mu = 3;
[T_A,X_A] = ExplicitEulerAdaptive(@VanDerPolJac,tspan,x0,h0,abstol,reltol,mu);

%Plotting
figure(2)
subplot(2,2,1)
plot(T_A,X_A(1,:),'b-','linewidth',2)
xlabel('t')
ylabel('$x_{1}$','Interpreter','Latex','Rotation',0,'FontSize',15)
title('$x_{1}(t)$','Interpreter','Latex','FontSize',18)

subplot(2,2,3)
plot(T_A,X_A(2,:),'b-','linewidth',2)
xlabel('t')
ylabel('$x_{2}$','Interpreter','Latex','Rotation',0,'FontSize',15)
title('$x_{2}(t)$','Interpreter','Latex','FontSize',18)

subplot(2,2,[2 4])
plot(X_A(1,:),X_A(2,:),'b-','linewidth',2)
xlabel('$x_{1}$','Interpreter','Latex','Rotation',0,'FontSize',15)
ylabel('$x_{2}$','Interpreter','Latex','Rotation',0,'FontSize',15)
title('Phase Plot','Interpreter','Latex','FontSize',18)

%% Van der pol Comparison with Matlab solvers. 
[TP_45,XP_45]=ode45(@VanDerPolJac,tspan,x0,[],mu);
[TP_23,XP_23]=ode23(@VanDerPolJac,tspan,x0,[],mu);
[TP_15s,XP_15s]=ode15s(@VanDerPolJac,tspan,x0,[],mu);

figure(11)
subplot(2,2,1)
plot(TP_45,XP_45(:,1))
hold on
plot(TP_23,XP_23(:,1))
plot(TP_15s,XP_15s(:,1))
legend('ode45','ode23','ode15s')
xlabel('t')
ylabel('x_1(t)')
hold off
subplot(2,2,3)
plot(TP_45,XP_45(:,2))
hold on
plot(TP_23,XP_23(:,2))
plot(TP_15s,XP_15s(:,2))
legend('ode45','ode23','ode15s')
xlabel('t')
ylabel('x_2(t)')
hold off
subplot(2,2,[2 4])
plot(XP_45(:,1),XP_45(:,2))
hold on
plot(XP_23(:,1),XP_23(:,2))
plot(XP_15s(:,1),XP_15s(:,2))
legend('ode45','ode23','ode15s')
xlabel('x_1(t)')
ylabel('x_2(t)')
hold off


%% Solve prey predator using Explicit Euler Fixed Tme Step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables;close all;clc

ta = 0;
tb = 25;

%Numerical Solution
N = 200;
x0 = [0.5 0.5];
a = 1;
b = 1;
p = [a b];
tspan = [ta tb];
[T,X] = ExplicitEuler(@preyPredator,tspan,N,x0,p);

%Plotting
figure(1)
subplot(2,2,1)
plot(T,X(1,:),'b-','linewidth',2)
xlabel('t')
ylabel('$x_{1}$','Interpreter','Latex','Rotation',0,'FontSize',15)
title('$x_{1}(t)$','Interpreter','Latex','FontSize',18)

subplot(2,2,3)
plot(T,X(2,:),'b-','linewidth',2)
xlabel('t')
ylabel('$x_{2}$','Interpreter','Latex','Rotation',0,'FontSize',15)
title('$x_{2}(t)$','Interpreter','Latex','FontSize',18)

subplot(2,2,[2 4])
plot(X(1,:),X(2,:),'b-','linewidth',2)
xlabel('$x_{1}$','Interpreter','Latex','Rotation',0,'FontSize',15)
ylabel('$x_{2}$','Interpreter','Latex','Rotation',0,'FontSize',15)
title('Phase Plot','Interpreter','Latex','FontSize',18)

%% Solve prey Predator using Explicit Euler Adaptive Time Step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Numerical Solution
tspan = [0 25];
x0 = [0.5 0.5];
h0 = 10^-2;
abstol = 10^(-2);
reltol = 10^(-2);
a = 1;
b = 1;
p = [a b];
[T_A,X_A] = ExplicitEulerAdaptive(@preyPredator,tspan,x0,h0,abstol,reltol,p);

%Plotting
figure(1)
subplot(2,2,1)
plot(T_A,X_A(1,:),'b-','linewidth',2)
xlabel('t')
ylabel('$x_{1}$','Interpreter','Latex','Rotation',0,'FontSize',15)
title('$x_{1}(t)$','Interpreter','Latex','FontSize',18)

subplot(2,2,3)
plot(T_A,X_A(2,:),'b-','linewidth',2)
xlabel('t')
ylabel('$x_{2}$','Interpreter','Latex','Rotation',0,'FontSize',15)
title('$x_{2}(t)$','Interpreter','Latex','FontSize',18)

subplot(2,2,[2 4])
plot(X_A(1,:),X_A(2,:),'b-','linewidth',2)
xlabel('$x_{1}$','Interpreter','Latex','Rotation',0,'FontSize',15)
ylabel('$x_{2}$','Interpreter','Latex','Rotation',0,'FontSize',15)
title('Phase Plot','Interpreter','Latex','FontSize',18)


%%
figure(1)
loglog([10^-8,10^-6,10^-4,10^-2],[111439,11182,1156,139],'o')
xlim([10^-9 10^-1])
xlabel('tol')
ylabel('steps')
%% Predator prey Comparison with Matlab solvers. 
[TP_45,XP_45]=ode45(@preyPredator,tspan,x0,[],p);
[TP_23,XP_23]=ode23(@preyPredator,tspan,x0,[],p);
[TP_15s,XP_15s]=ode15s(@preyPredator,tspan,x0,[],p);

figure(11)
subplot(2,2,1)
plot(TP_45,XP_45(:,1))
hold on
plot(TP_23,XP_23(:,1))
plot(TP_15s,XP_15s(:,1))
legend('ode45','ode23','ode15s')
hold off
subplot(2,2,3)
plot(TP_45,XP_45(:,2))
hold on
plot(TP_23,XP_23(:,2))
plot(TP_15s,XP_15s(:,2))
legend('ode45','ode23','ode15s')
hold off
subplot(2,2,[2 4])
plot(XP_45(:,1),XP_45(:,2))
hold on
plot(XP_23(:,1),XP_23(:,2))
plot(XP_15s(:,1),XP_15s(:,2))
legend('ode45','ode23','ode15s')
hold off