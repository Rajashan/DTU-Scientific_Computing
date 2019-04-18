%% For testing stability of solver on stiff stiff problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables;close all;clc

ta = 0;
tb = 3;

%Numerical Solution
N = 14;
x0 = 1;
lambda = -10;
tspan = [ta tb];
[T,X_I] = ImplicitEuler(@TestEquationJac,tspan,N,x0,lambda);
[~,X_E] = ExplicitEuler(@TestEquationJac,tspan,N,x0,lambda);

%Analytical Solution
X_A = x0*exp(lambda*T);

figure(1)
plot(T,X_A)
hold on
plot(T,X_E,'-o')
plot(T,X_I,'-o')
hold off
title('dt=0.2')
legend('Analytical','Explicit','Implicit')
xlabel('t')
ylabel('x(t)')