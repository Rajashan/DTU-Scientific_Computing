clear;
close all;
%% Van der Pol, constants and std Wiener

T  =25;
N = 300;
Ns = 10;
seed = 1000;
Nw = 1;


%% Wiener Process 
[W,TW,~] = StdWienerProcess(T,N,Nw,Ns,seed);


%% Explicit Explicit Van der Pol
mu = 10;
sigma = 1;
p = [mu sigma];
x0 = [0.5 0.5];
X_V1 = zeros(length(x0),N+1,Ns);
for i=1:Ns
    X_V1(:,:,i) = SDE_ExplicitExplicitEuler(...
        @VanderpolDrift,@VanderpolDiffusion1,...
        TW,x0,W(:,:,i),p);
end

X_A1(:,:) = SDE_ExplicitExplicitEuler(...
        @VanderpolDrift,@VanderpolDiffusion1,...
        TW,x0,W(:,:),[mu 0.0]);
   
[xmean_V1,s_V1,xmeanp2s_V1,xmeanm2s_V1]=ScalarSampleMeanStdVar(X_V1);
X_V2 = zeros(length(x0),N+1,Ns);

for i=1:Ns
    X_V2(:,:,i) = SDE_ExplicitExplicitEuler(...
        @VanderpolDrift,@VanderpolDiffusion2,...
        TW,x0,W(:,:,i),p);
end

X_A2(:,:) = SDE_ExplicitExplicitEuler(...
        @VanderpolDrift,@VanderpolDiffusion2,...
        TW,x0,W(:,:),[mu 0.0]);
figure(1)
title('Explicit explicit, diffusion 1')
subplot(2,2,1)
plot(TW,squeeze(X_V1(1,:,:)))
hold on
plot(TW,squeeze(X_A1(1,:,:)),'linewidth',2,'color','k')
hold off
xlabel('t')
ylabel('x_1')
subplot(2,2,2)
plot(TW,squeeze(X_V1(2,:,:)))
hold on
plot(TW,squeeze(X_A1(2,:,:)),'linewidth',2,'color','k')
hold off
xlabel('t')
ylabel('x_2')
subplot(2,2,[3 4])
plot(squeeze(X_V1(1,:,:)),squeeze(X_V1(2,:,:)))
hold on
plot(squeeze(X_A1(1,:,:)),squeeze(X_A1(2,:,:)),'linewidth',2,'color','k')
hold off
xlabel('x_1')
ylabel('x_2')

figure(2)
subplot(2,2,1)
plot(TW,squeeze(X_V2(1,:,:)))
hold on
plot(TW,squeeze(X_A2(1,:,:)),'linewidth',2,'color','k')
hold off
xlabel('t')
ylabel('x_1')
subplot(2,2,2)
plot(TW,squeeze(X_V2(2,:,:)))
hold on
plot(TW,squeeze(X_A2(2,:,:)),'linewidth',2,'color','k')
hold off
xlabel('t')
ylabel('x_2')
subplot(2,2,[3 4])
plot(squeeze(X_V2(1,:,:)),squeeze(X_V2(2,:,:)))
hold on
plot(squeeze(X_A2(1,:,:)),squeeze(X_A2(2,:,:)),'linewidth',2,'color','k')
hold off
xlabel('x_1')
ylabel('x_2')
%% Explicit Implicit Van der Pol
X2_V1 = zeros(length(x0),N+1,Ns);
for i=1:Ns
    X2_V1(:,:,i) = SDE_ImplicitExplicitEuler(...
        @VanderpolDrift,@VanderpolDiffusion1,...
        TW,x0,W(:,:,i),p);
end

X_I1(:,:) = SDE_ImplicitExplicitEuler(...
        @VanderpolDrift,@VanderpolDiffusion1,...
        TW,x0,W(:,:),[mu 0.0]);

X2_V2 = zeros(length(x0),N+1,Ns);
for i=1:Ns
    X2_V2(:,:,i) = SDE_ImplicitExplicitEuler(...
        @VanderpolDrift,@VanderpolDiffusion2,...
        TW,x0,W(:,:,i),p);
end

X_I2(:,:) = SDE_ImplicitExplicitEuler(...
        @VanderpolDrift,@VanderpolDiffusion2,...
        TW,x0,W(:,:),[mu 0.0]);
figure(1)
subplot(2,2,1)
plot(TW,squeeze(X2_V1(1,:,:)))
hold on
plot(TW,squeeze(X_I1(1,:,:)),'linewidth',2,'color','k')
hold off
xlabel('t')
ylabel('x_1')
subplot(2,2,2)
plot(TW,squeeze(X2_V1(2,:,:)))
hold on
plot(TW,squeeze(X_I1(2,:,:)),'linewidth',2,'color','k')
hold off
xlabel('t')
ylabel('x_2')
subplot(2,2,[3 4])
plot(squeeze(X2_V1(1,:,:)),squeeze(X2_V1(2,:,:)))
hold on
plot(squeeze(X_I1(1,:,:)),squeeze(X_I1(2,:,:)),'linewidth',2,'color','k')
hold off
xlabel('x_1')
ylabel('x_2')

figure(2)
subplot(2,2,1)
plot(TW,squeeze(X2_V2(1,:,:)))
hold on
plot(TW,squeeze(X_I2(1,:,:)),'linewidth',2,'color','k')
hold off
xlabel('t')
ylabel('x_1')
subplot(2,2,2)
plot(TW,squeeze(X2_V2(2,:,:)))
hold on
plot(TW,squeeze(X_I2(2,:,:)),'linewidth',2,'color','k')
hold off
xlabel('t')
ylabel('x_2')
subplot(2,2,[3 4])
plot(squeeze(X2_V2(1,:,:)),squeeze(X2_V2(2,:,:)))
hold on
plot(squeeze(X_I2(1,:,:)),squeeze(X_I2(2,:,:)),'linewidth',2,'color','k')
hold off
xlabel('x_1')
ylabel('x_2')

%% Exact sol


[xmean3,s3,xmeanp2s3,xmeanm2s3]=ScalarSampleMeanStdVar(exact);
figure(10)
plot(TW,xmean)
hold on
plot(TW,xmeanp2s)
plot(TW,xmeanm2s)


figure(20)
plot(TW,xmean2)
hold on
plot(TW,xmeanp2s2)
plot(TW,xmeanm2s2)
hold off

figure(30)
plot(TW,xmean3)
hold on
plot(TW,xmeanp2s3)
plot(TW,xmeanm2s3)
hold off

figure(40)
histogram(X(1,N+1,:))
figure(50)
histogram(X2(1,N+1,:))
figure(60)
histogram(exact(1,N+1,:))

