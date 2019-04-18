clear all;
close all;

theta= linspace(0,2*pi,1000); % The boundary interval
x = exp(i*theta); % The boundary interval in x
z = (x.^2-1)./(1/3.*x.^2+4/3.*x+1/3); % Combined polynomial 
% Plot boundary of area of absolute stability
figure()
plot(real(z),imag(z),'b-','linewidth',3)
title('Region of absolute stability')
lgd = legend('Region boundary');
xlabel('Real axis')
ylabel('Imaginary axis')
xlim([-1,1]);
set(gca,'fontsize',20)
lgd.FontSize = 18;