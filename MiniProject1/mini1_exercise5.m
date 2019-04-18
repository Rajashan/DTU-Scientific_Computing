clear;
close all;
%% Geometric Brownian, constants and std Wiener

T  =10;
N = 100;
Ns = 10;
seed = 1000;
Nw = 1;
x0 = 1;
lambda = 0.15;
sigma = 0.15;
p = [lambda sigma];

%% Wiener Process
[W,TW,~] = StdWienerProcess(T,N,Nw,Ns,seed);


%% Explicit Explicit
X = zeros(length(x0),N+1,Ns);
for i=1:Ns
    X(:,:,i) = SDE_ExplicitExplicitEuler(...
        @GeometricBrownianDrift,@GeometricBrownianDiffusion,...
        TW,x0,W(:,:,i),p);
end

[xmean,s,xmeanp2s,xmeanm2s]=ScalarSampleMeanStdVar(X);


%% Explicit Implicit
X2 = zeros(length(x0),N+1,Ns);
for i=1:Ns
    X2(:,:,i) = SDE_ImplicitExplicitEuler(...
        @GeometricBrownianDrift,@GeometricBrownianDiffusion,...
        TW,x0,W(:,:,i),p);
end

[xmean2,s2,xmeanp2s2,xmeanm2s2]=ScalarSampleMeanStdVar(X2);



%% Exact sol
exact = x0*exp((lambda-(1/2)*sigma^2)*TW+sigma*W);


[xmean3,s3,xmeanp2s3,xmeanm2s3]=ScalarSampleMeanStdVar(exact);


%% Plots + histogram of ex, im, exact

figure(1)
subplot(1,3,1)
plot(TW,squeeze(X(1,:,:)))
hold on
%plot(TW,xmean,'linewidth',2,'color','k')
%plot(TW,xmeanp2s,'linewidth',2,'color','b')
%plot(TW,xmeanm2s,'linewidth',2,'color','r')
hold off
title('Explicit-Explicit Euler-Maruyama, 10000 realizations')
xlabel('t')
ylabel('x(t)')
legend('mean','+2std','-2std')
subplot(1,3,2)
plot(TW,squeeze(X2(1,:,:)))
hold on
%plot(TW,xmean2,'linewidth',2,'color','k')
%plot(TW,xmeanp2s2,'linewidth',2,'color','b')
%plot(TW,xmeanm2s2,'linewidth',2,'color','r')
hold off
title('Implicit-Explicit Euler-Maruyama, 10000 realizations')
xlabel('t')
ylabel('x(t)')
legend('mean','+2std','-2std')
subplot(1,3,3)
plot(TW,squeeze(exact(1,:,:)))
hold on
%plot(TW,xmean3,'linewidth',2,'color','k')
%plot(TW,xmeanp2s3,'linewidth',2,'color','b')
%plot(TW,xmeanm2s3,'linewidth',2,'color','r')
hold off
title('Exact solution, 10000 realizations')
xlabel('t')
ylabel('x(t)')
legend('mean','+2std','-2std')
%%
figure(2)
subplot(1,3,1)
histogram(X(1,N+1,:),20)
title('Explicit-Explicit Euler-Maruyama, 10000 realizations')
xlabel('x(t_f)')
subplot(1,3,2)
histogram(X2(1,N+1,:),20)
title('Implicit-Explicit Euler-Maruyama, 10000 realizations')
xlabel('x(t_f)')
subplot(1,3,3)
histogram(exact(1,N+1,:),20)
title('Exact solution, 10000 realizations')
xlabel('x(t_f)')
