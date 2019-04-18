clear all;
close all;
%% 9-point Laplacian, fourth order accurate
a = 0; b = 1; c = 0; d = 1; % Define region
m = 100; % number of grid points
n = m+2; % number of grid points with boundaries
h = 1/(m+1); % grid point size
[X,Y] = meshgrid(a+h:h:b-h,c+h:h:d-h); % Define grid without boundaries
[bcX,bcY] = meshgrid(a:h:b,c:h:d); % Define grid with boundaries 

u0 = @(x,y) real(sin(2*pi*(x-y).^(2.5)));% real(sin(2*pi*(x-y).^(2.5))) + imag(sin(2*pi*(x-y).^(2.5))); % Define exact solution
uTrue = u0(bcX,bcY) + u0(bcY,bcX); % Define true solution
border = [ones(1,n);[ones(n-2,1),zeros(n-2,n-2),ones(n-2,1)];ones(1,n)]; % Define border
bc = (u0(bcX,bcY) + u0(bcY,bcX)).*border; % Define baoundary at border
bcEdge = bc; % Define the edge as the boundary
bcEdge(bc == 0) = NaN; % Remove every other point than the boundary

fun = @(x,y) 47.12388980*(abs(x-y)).^(1/2).*sign(x-y).^2.*cos(2*pi*abs(x-y).^(2.5)) +... 
           0*31.41592654*(abs(x-y)).^(3/2).*sign(x-y)   .*cos(2*pi*abs(x-y).^(2.5)) -...
             493.4802202*(abs(x-y)).^(3/1).*sign(x-y).^2.*sin(2*pi*abs(x-y).^(2.5));

fOld = zeros(m^2,1);
for ii = 1:m
    fOld(1+(ii-1)*m:ii*m,1) = fun(X(ii,:),Y(ii,1));
end
LapRhs = LaplacianOfRhs(m,h,fun,bcX,bcY);
f = fOld + h^2/12*LapRhs;

A = poisson9(m);
F = form_rhs9(m,f,bc);

U = A\F; % Calculate the solution U (as a vector)
UU = zeros(m,m); 
for ii = 1:m % Reshape U as a m x m matrix
    UU(ii,:) = U(1+(ii-1)*m:ii*m,1);
end

subplot(1,2,1)
hold on;
s = surf(bcX,bcY,uTrue);
s.EdgeColor = 'none';
surf(bcX,bcY,bcEdge,'EdgeColor','r');
title('True solution for u_{2,exact}')
xlabel('x-axis');
ylabel('y-axis');
set(gca,'fontsize',20);
subplot(1,2,2)
hold on;
s = surf(X,Y,UU);
s.EdgeColor = 'none';
surf(bcX,bcY,bcEdge,'EdgeColor','r');
title('Solution generated with 9-point Laplacian for u_{2,exact}') 
xlabel('x-axis');
ylabel('y-axis');
zlim([ -1.1 1.1])
set(gca,'fontsize',20);

%% Convergence for the 9-point Laplacian, fourth order Accuracy

for kk = 3:10
    h = 1/2^(kk); % grid point size
    hx(kk-2,1) = h;
    m = -1+1/h; % number of grid points
    n = m+2; % number of grid points with boundaries
    [X,Y] = meshgrid(a+h:h:b-h,c+h:h:d-h); % Define grid without boundaries
    [bcX,bcY] = meshgrid(a:h:b,c:h:d); % Define grid with boundaries 
    
    uTrue = u0(bcX,bcY) + u0(bcY,bcX); % Define true solution
    border = [ones(1,n);[ones(n-2,1),zeros(n-2,n-2),ones(n-2,1)];ones(1,n)]; % Define border
    bc = (u0(bcX,bcY) + u0(bcY,bcX)).*border; % Define baoundary at border
    bcEdge = bc; % Define the edge as the boundary
    bcEdge(bc == 0) = NaN; % Remove every other point than the boundary

    fun = @(x,y) 47.12388980*(abs(x-y)).^(1/2).*sign(x-y).^2.*cos(2*pi*abs(x-y).^(2.5)) +... 
               0*31.41592654*(abs(x-y)).^(3/2).*sign(x-y)   .*cos(2*pi*abs(x-y).^(2.5)) -...
                 493.4802202*(abs(x-y)).^(3/1).*sign(x-y).^2.*sin(2*pi*abs(x-y).^(2.5));
    
    fOld = zeros(m^2,1);
    for ii = 1:m
        fOld(1+(ii-1)*m:ii*m,1) = fun(X(ii,:),Y(ii,1));
    end
    LapRhs = LaplacianOfRhs(m,h,fun,bcX,bcY);
    f = fOld + h^2/12*LapRhs;

    A = poisson9(m);
    F = form_rhs9(m,f,bc);

    U = A\F; % Calculate the solution U (as a vector)
    UU = zeros(m,m); 
    for ii = 1:m % Reshape U as a m x m matrix
        UU(ii,:) = U(1+(ii-1)*m:ii*m,1);
    end
    
    GTE(kk-2,1) = max(max(abs(UU-uTrue(2:m+1,2:m+1)))); % Infinity norm 
end

%% Plot rate of convergence for 9-point Laplacian
figure;
hold on;
plot(hx,hx.^(3/2),'k-','linewidth',2)
plot(hx,GTE,'r--*','linewidth',2)
title('Rate of convergence for 9-point Laplacian with u_{2,exact} as solution')
ylabel('Global truncation error')
xlabel('stepsize, h')
legend('Order of 3/2 (h^{3/2})','Actual convergence','location','north')
grid on;
set(gca,'fontsize',20, 'XScale', 'log', 'YScale', 'log')
xlim([min(hx)*0.9 max(hx)*1.1])
