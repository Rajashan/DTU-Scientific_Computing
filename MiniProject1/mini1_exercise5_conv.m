%% Exercise 9
clear variables; close all; clc
%% Wiener Process
close all;
Ns = 10.^5;    %Realizations or Number of Particles
T = 10;     %Time 0 to T
N = 1000;   %Number of Steps
seed = 100;   %seed
[W,Tw,dW,S] = WienerProcess(Ns,N,T,seed);

%WARNING PLOTTING MAY CRASH MATLAB OR TAKE VERY LONG TIME REDUCE Ns TO
%APPROPRIATE NUMBER OF REALIZATIONS INSTEAD OF 10^5
%{ 
subplot(2,1,1)
plot(Tw,W,'linewidth',2);
title('$W(t)$','Interpreter','Latex','Fontsize',18)

subplot(2,1,2)
plot(Tw(1:end-1),dW,'linewidth',2);
title('$dW(t)$','Interpreter','Latex','Fontsize',18)
%}

figure()
plot(Tw,W(1:10,:),'linewidth',2);
title('$W(t)$ for ten first realizations','Interpreter','Latex','Fontsize',18)
hold on
plot(Tw,zeros(1,length(Tw)),'linewidth',3,'Color','red')
plot(Tw,[2*sqrt(Tw) ; -2*sqrt(Tw)],'linewidth',3,'Color','red')
plot(Tw,S,'linewidth',3,'Color','black')

%% Example 1: Constant Coefficient - SDE_ExplicitExplicitEuler
close all
T = 10;
N = 1000;
Ns = 10^5;
seed = 100;

lambda = 1; 
sigma = 1;
x0 = 0;
fun = @(t,x) lambda;
gun = @(t,x) sigma;

[W,Tw,Dw,S] = WienerProcess(Ns,N,T,seed);
X = zeros(size(W));
for i = 1:10
X(i,:) = SDE_ExplicitExplicitEuler(fun,gun,Tw,x0,W(i,:));
end
for i = 1:3
X(i+10,:) = SDE_ExplicitExplicitEuler(fun,gun,Tw,x0,S(i,:));
end

figure(1)
title({'Constant Coefficients','$dx(t) = \lambda dt + \sigma d\omega(t)$'},'Interpreter','Latex','Fontsize',18)
hold on
plot(Tw,X(1:10,:),'linewidth',3)
plot(Tw,X(11:13,:),'Linewidth',3,'Color','black')
ylim([-4 16])



%% Example 2: Geometric Brownian Motion - - SDE_ExplicitExplicitEuler
clear variables; close all
T = 10;
N = 1000;
Ns = 10^5;      %Increase: Black lines to convergence towards mean and std.
seed = 100;

lambda = 0.15;
sigma = 0.15;
x0 = 1;
fun = @(t,x) lambda*x;
gun = @(t,x) sigma*x;

[W,Tw,Dw,S] = WienerProcess(Ns,N,T,seed);
for i = 1:10
X(i,:) = SDE_ExplicitExplicitEuler(fun,gun,Tw,x0,W(i,:));
end
for i = 1:3
X(i+10,:) = SDE_ExplicitExplicitEuler(fun,gun,Tw,x0,S(i,:));
end

figure(1)
title({'Geometric Brownian','$dx(t) = \lambda x(t) dt + \sigma x(t) d\omega(t)$'},'Interpreter','Latex','Fontsize',18)
hold on
plot(Tw,X(1:10,:),'linewidth',3)
plot(Tw,X(11:13,:),'Linewidth',3,'Color','black')
ylim([0 14])

%% Calculating the Order of Convergence
T = 10;
N = 10.^(1:7);
Ns = 10^0;
seed = 100;

lambda = 0.15;
sigma = 0.15;
p = [lambda,sigma];
x0 = 1;

for ii = 1:length(N)
    [W,Tw,Dw,S] = WienerProcess(Ns,N(ii),T,seed);
    X_Exp = SDE_ExplicitExplicitEuler(@BrownianMotionJac,@BrownianMotionSDE,Tw,x0,W,p);
    X_Imp = SDE_ImplicitExplicitEuler(@BrownianMotionJac,@BrownianMotionSDE,Tw,x0,W,p);
    Y = x0*exp( (lambda - 0.5*sigma^2)*Tw + sigma*W);
    E1(ii) = mean(X_Exp)-mean(Y);
    E2(ii) = mean(X_Imp)-mean(Y);
end

%%
h = T./N;
%% Demonstrates that the order of the methods is p = 1/2
figure(1)
plot(h,E1,'LineStyle','--','LineWidth',3,'Color','Black')
hold on
plot(h,E2,'LineStyle','--','LineWidth',3,'Color','Red')
plot(h,h,'Linewidth',3)
set(gca,'XScale','log','YScale','log')
title('Convergence Plot for SDE, 10.000 realizations','Interpreter','Latex','FontSize',18)
xlabel('Time Step h')
ylabel('Error')
lg = legend('Stochastic ExplicitExplicitEuler Global Truncation Error','Stochastic ImplicitExplicitEuler Global Truncation Error','Reference Line Order $N^{1/2}$');
set(lg,'Interpreter','Latex','FontSize',15,'location','southeast')
ylim([10^-7 10^2])
hold off

%% Example 3: Langevin - - SDE_ExplicitExplicitEuler
close all
T = 10;
N = 1000;
Ns = 10^5;      %Increase: Black lines to convergence towards mean and std.
seed = 100;

lambda = -0.5;
sigma = 1;
x0 = 10;
fun = @(t,x) lambda*x;
gun = @(t,x) sigma;

[W,Tw,~,S] = WienerProcess(Ns,N,T,seed);
X = zeros(size(W));
for i = 1:10
X(i,:) = SDE_ExplicitExplicitEuler(fun,gun,Tw,x0,W(i,:));
end
for i = 1:3
X(i+10,:) = SDE_ExplicitExplicitEuler(fun,gun,Tw,x0,S(i,:));
end

figure(1)
title({'Langevin','$dx(t) = \lambda x(t) dt + \sigma d\omega(t)$'},'Interpreter','Latex','Fontsize',18)
hold on
plot(Tw,X(1:10,:),'linewidth',3)
plot(Tw,X(11:13,:),'Linewidth',3,'Color','black')
ylim([-4 12])

%% Multidimensional Wiener Process
close all; clear variables;

T = 10;
N = 1000;
nW = 2;
Ns = 10^4;
seed = 100;

[W,Tw,dW,S] = MultiWienerProcess(T,N,nW,Ns,seed);

figure()
for ii = 1:10
plot3(W(1,:,ii),W(2,:,ii),Tw,'linewidth',2);
hold on
end
xlabel('T')
ylabel('W')
title('$W(t)$ for ten first realizations','Interpreter','Latex','Fontsize',18)
hold on
plot3(S(1,:),S(2,:),Tw,'linewidth',4,'Color','Black');
plot3(S(3,:),S(4,:),Tw,'linewidth',4,'Color','Black');
plot3(S(5,:),S(6,:),Tw,'linewidth',4,'Color','Black');
%% Van Der Pol - SDE_ExplicitExplicitEuler
clear variables; close all;

T = 15;
N = 1000;
Ns = 10^4;
Nw = 2;
seed = 100;

[W,Tw,dW,S] = MultiWienerProcess(T,N,Nw,Ns,seed);

mu = 3;
sigma = 0.2;
p = [mu, sigma];
x0 = [0.5 0.5];

for i = 1:10
    X(:,:,i) = SDE_ExplicitExplicitEuler(@VanDerPolJac,@VanDerPolSDE,Tw,x0,W(:,:,i),p);
end
%Xr is reference solution to pure ODE sigma = 0 (without stochastic term)
Xr = SDE_ExplicitExplicitEuler(@VanDerPolJac,@VanDerPolSDE,Tw,x0,W(:,:,i),[mu,0]);
Xs = SDE_ExplicitExplicitEuler(@VanDerPolJac,@VanDerPolSDE,Tw,x0,S(1:2,:),p);
Xsp2 = SDE_ExplicitExplicitEuler(@VanDerPolJac,@VanDerPolSDE,Tw,x0,S(3:4,:),p);
Xsm2 = SDE_ExplicitExplicitEuler(@VanDerPolJac,@VanDerPolSDE,Tw,x0,S(5:6,:),p);

figure(1)
subplot(2,2,1)
for ii = 1:10
    plot(Tw,X(1,:,ii),'Linewidth',2)
    hold on
end
plot(Tw,Xs(1,:),'Linewidth',4,'Color','Black')
%plot(Tw,Xsp2(1,:),'Linewidth',4,'Color','Black')
%plot(Tw,Xsm2(1,:),'Linewidth',4,'Color','Black')
plot(Tw,Xr(1,:),'Linewidth',4,'Color','Red','Linestyle','--')
ylim([-3 3])
xlabel('T')
ylabel('$x_{1}(t)$','Interpreter','Latex','FontSize',15)
title('VanDerPol solution $x_{1}(t)$','Interpreter','Latex','FontSize',15)

figure(1)
subplot(2,2,2)
for ii = 1:10
    plot(Tw,X(2,:,ii),'Linewidth',2)
    hold on
end
plot(Tw,Xs(2,:),'Linewidth',4,'Color','Black')
plot(Tw,Xsp2(2,:),'Linewidth',4,'Color','Black')
plot(Tw,Xsm2(2,:),'Linewidth',4,'Color','Black')
plot(Tw,Xr(2,:),'Linewidth',4,'Color','Red','Linestyle','--')
ylim([-7 7])
xlabel('T')
ylabel('$x_{2}(t)$','Interpreter','Latex','FontSize',15)
title('VanDerPol solution $x_{2}(t)$','Interpreter','Latex','FontSize',15)

%figure(3)
subplot(2,2,3:4)
for ii = 1:10
    plot(X(1,:,ii),X(2,:,ii),'Linewidth',2)
    hold on
end
plot(Xs(1,:),Xs(2,:),'Linewidth',3,'Color','Black')
plot(Xsp2(1,:),Xsp2(2,:),'Linewidth',3,'Color','Black')
plot(Xsm2(1,:),Xsm2(2,:),'Linewidth',3,'Color','Black')
ylim([-8 8])
xlim([-2.5 2.5])
xlabel('$x_{1}(t)$','Interpreter','Latex','FontSize',15)
ylabel('$x_{2}(t)$','Interpreter','Latex','FontSize',15)
title('VanDerPol $x_{1}(t)$/$x_{2}(t)$','Interpreter','Latex','FontSize',15)
set(gcf,'units','points','position',[150,0,1000,1000])

%% Calculating Order of Method
clear variables; close all;

T = 15;
N = 5.^(3:9);
Ns = 1;
Nw = 2;
seed = 100;

mu = 3;
sigma = 0.2;
p = [mu, sigma];
x0 = [0.5 0.5];

%Takes 40 seconds

for ii = 1:length(N)
    tic
    [W,Tw,dW,S] = MultiWienerProcess(T,N(ii),Nw,Ns,seed);
    X{ii} = SDE_ExplicitExplicitEuler(@VanDerPolJac,@VanDerPolSDE,Tw,x0,S(1:2,:),p);
    Xr{ii} = SDE_ExplicitExplicitEuler(@VanDerPolJac,@VanDerPolSDE,Tw,x0,W(:,:),[mu,0]);
    toc
end

%Expected Value

for ii = 1:length(X)
   E(:,ii) = mean(abs(Xr{ii} - X{ii}),2);
end

figure(1)
plot(N,E(1,:))
set(gca,'XScale','log','YScale','log')

%figure(2)

%% Brownian Motion - SDE_ExplicitExplicitEuler

