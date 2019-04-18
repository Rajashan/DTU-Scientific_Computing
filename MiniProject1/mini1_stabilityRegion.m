clear;
close all;
% Absolute stability region plot. 

theta= linspace(0,2*pi,1000); % The boundary interval
x = exp(1i*theta); % The boundary interval in x
z = 1+x; % Combined polynomial 
% Plot boundary of area of absolute stability
figure(1)
plot(real(z),imag(z),'k-','linewidth',2)
title('Region of absolute stability, Implicit Euler')
lgd = legend('Region boundary');
xlabel('Real axis')
ylabel('Imaginary axis')
xlim([-4,4]);
ylim([-4,4]);
set(gca,'fontsize',20)
lgd.FontSize = 18;
xticks([-4 -2 0 2 4])
yticks([-4 -1  1 4])
grid on