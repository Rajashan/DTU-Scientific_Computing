clear all
close all;
%% 5-point Laplacian

% Define region
a = 0;
b = 1;
c = 0;
d = 1;

m = 100; % number of grid points
n = m+2; % number of grid points with boundaries
h = 1/(m+1); % grid point size
[X,Y] = meshgrid(a+h:h:b-h,c+h:h:d-h); % Define grid without boundaries
[bcX,bcY] = meshgrid(a:h:b,c:h:d); % Define grid with boundaries 

u0 = @(x,y) sin(4*pi*(x.^2 + y.^2)) + cos(4*pi*x.*y); % Define exact solution
uTrue = u0(bcX,bcY); % Define true solution
border = [ones(1,n);[ones(n-2,1),zeros(n-2,n-2),ones(n-2,1)];ones(1,n)]; % Define border
bc = u0(bcX,bcY).*border; % Define baoundary at border
bcEdge = bc; % Define the edge as the boundary
bcEdge(bc == 0) = NaN; % Remove every other point than the boundary

% figure; % Plot true solution and boundary
% subplot(1,2,1)
% hold on;
% s = surf(bcX,bcY,uTrue);
% s.EdgeColor = 'none';
% surf(bcX,bcY,bcEdge,'EdgeColor','r');
% title('True solution')
% subplot(1,2,2)
% surf(bcX,bcY,bc)
% title('Boundary condition')

A = poisson5(m); % Calle function to create A spare matrix
% full(A)
% Define the function f as u''(x,y)
fun = @(x,y) 16*pi.*(cos(4*pi*(x.^2+y.^2))-pi.*(x.^2+y.^2).*(4.*sin(4*pi*(x.^2+y.^2))+cos(4*pi*x.*y)));
f = zeros(m^2,1); % Initialiez f as a m^2 x 1 vector
for ii = 1:m % Create the value for the f-vector
    f(1+(ii-1)*m:ii*m,1) = fun(X(ii,:),Y(ii,1));
end

F = form_rhs5(m ,f,bc); % Form the right hand side

U = A\F; % Calculate the solution U (as a vector)
UU = zeros(m,m); 
for ii = 1:m % Reshape U as a m x m matrix
    UU(ii,:) = U(1+(ii-1)*m:ii*m,1);
end
% Plot solution
figure;
subplot(1,2,1)
hold on;
s = surf(bcX,bcY,uTrue);
s.EdgeColor = 'none';
surf(bcX,bcY,bcEdge,'EdgeColor','r');
title('True solution')
subplot(1,2,2)
hold on;
s = surf(X,Y,UU);
% s = surf(X,Y,UU-uTrue(1:m,1:m));
s.EdgeColor = 'none';
surf(bcX,bcY,bcEdge,'EdgeColor','r');
title('Approximated solution')

%% 9-point Laplacian

a = 0; b = 1; c = 0; d = 1; % Define region
m = 100; % number of grid points
n = m+2; % number of grid points with boundaries
h = 1/(m+1); % grid point size
[X,Y] = meshgrid(a+h:h:b-h,c+h:h:d-h); % Define grid without boundaries
[bcX,bcY] = meshgrid(a:h:b,c:h:d); % Define grid with boundaries 

u0 = @(x,y) sin(4*pi*(x.^2 + y.^2)) + cos(4*pi*x.*y); % Define exact solution
uTrue = u0(bcX,bcY); % Define true solution
border = [ones(1,n);[ones(n-2,1),zeros(n-2,n-2),ones(n-2,1)];ones(1,n)]; % Define border
bc = u0(bcX,bcY).*border; % Define baoundary at border
bcEdge = bc; % Define the edge as the boundary
bcEdge(bc == 0) = NaN; % Remove every other point than the boundary

A = poisson9(m);

fun = @(x,y) 16*pi.*(cos(4*pi*(x.^2+y.^2))-pi.*(x.^2+y.^2).*(4.*sin(4*pi*(x.^2+y.^2))+cos(4*pi*x*y)));
f = zeros(m^2,1);
for ii = 1:m
    f(1+(ii-1)*m:ii*m,1) = fun(X(ii,:),Y(ii,1));
end

F = form_rhs9(m,f,bc);

U = A\F; % Calculate the solution U (as a vector)
UU = zeros(m,m); 
for ii = 1:m % Reshape U as a m x m matrix
    UU(ii,:) = U(1+(ii-1)*m:ii*m,1);
end

% Plot solution
figure;
subplot(1,2,1)
hold on;
s = surf(bcX,bcY,uTrue);
s.EdgeColor = 'none';
surf(bcX,bcY,bcEdge,'EdgeColor','r');
title('True solution')
subplot(1,2,2)
hold on;
s = surf(X,Y,UU);
s.EdgeColor = 'none';
surf(bcX,bcY,bcEdge,'EdgeColor','r');
title('Approximated solution')

%% Convergence for standard 9-point Laplacian

a = 0; b = 1; c = 0; d = 1;
u0 = @(x,y) sin(4*pi*(x.^2 + y.^2)) + cos(4*pi*x.*y); % Define exact solution
fun = @(x,y) 16*pi.*(cos(4*pi*(x.^2+y.^2))-pi.*(x.^2+y.^2).*(4.*sin(4*pi*(x.^2+y.^2))+cos(4*pi*x.*y)));
                       
for kk = 3:10
    
    h = 1/2^(kk); % grid point size
    hx1(kk-2,1) = h;
    m = -1+1/h; % number of grid points
    n = m+2; % number of grid points with boundaries
    
    [X,Y] = meshgrid(a+h:h:b-h,c+h:h:d-h); % Define grid without boundaries
    [bcX,bcY] = meshgrid(a:h:b,c:h:d); % Define grid with boundaries 

    uTrue = u0(bcX,bcY); % Define true solution
    border = [ones(1,n);[ones(n-2,1),zeros(n-2,n-2),ones(n-2,1)];ones(1,n)]; % Define border
    bc = u0(bcX,bcY).*border; % Define baoundary at border
    
    A = poisson9(m);

    f = zeros(m^2,1);
    for ii = 1:m
        f(1+(ii-1)*m:ii*m,1) = fun(X(ii,:),Y(ii,1));
    end

    F = form_rhs9(m,f,bc);

    U = A\F; % Calculate the solution U (as a vector)
    UU = zeros(m,m);
    for ii = 1:m % Reshape U as a m x m matrix
        UU(ii,:) = U(1+(ii-1)*m:ii*m,1);
    end
    GTE1(kk-2,1) = max(max(abs(UU-uTrue(2:m+1,2:m+1)))); % Infinity norm
end

%% Plot rate of convergence for 9-point Laplacian
figure;
hold on;
plot(hx1,hx1.^2,'k-','linewidth',2)
plot(hx1,GTE1,'r--*','linewidth',2)
title('Rate of convergence of 9-point Laplacian')
ylabel('Global truncation error')
xlabel('stepsize, h')
legend('Expected Second order','Actual convergence','location','north')
grid on;
set(gca,'fontsize',20, 'XScale', 'log', 'YScale', 'log')

%% 9-point Laplacian, fourth order accurate

a = 0; b = 1; c = 0; d = 1;
m = 100; % number of grid points
n = m+2; % number of grid points with boundaries
h = 1/(m+1); % grid point size
[X,Y] = meshgrid(a+h:h:b-h,c+h:h:d-h); % Define grid without boundaries
[bcX,bcY] = meshgrid(a:h:b,c:h:d); % Define grid with boundaries 
u0 = @(x,y) sin(4*pi*(x.^2 + y.^2)) + cos(4*pi*x.*y); % Define exact solution
uTrue = u0(bcX,bcY); % Define true solution
border = [ones(1,n);[ones(n-2,1),zeros(n-2,n-2),ones(n-2,1)];ones(1,n)]; % Define border

fun = @(x,y) 16*pi.*(cos(4*pi*(x.^2+y.^2))-pi.*(x.^2+y.^2).*(4.*sin(4*pi*(x.^2+y.^2))+cos(4*pi*x.*y)));

fOld = zeros(m^2,1);
for ii = 1:m
    fOld(1+(ii-1)*m:ii*m,1) = fun(X(ii,:),Y(ii,1));
end
LapRhs = LaplacianOfRhs(m,h,fun,bcX,bcY);
f = fOld + h^2/12*LapRhs;

bc = u0(bcX,bcY).*border; % Define baoundary at border
bcEdge = bc; % Define the edge as the boundary
bcEdge(bc == 0) = NaN; % Remove every other point than the bounda

A9 = poisson9(m);
F = form_rhs9(m,f,bc);

U = A9\F; % Calculate the solution U (as a vector)
UU = zeros(m,m); 
for ii = 1:m % Reshape U as a m x m matrix
    UU(ii,:) = U(1+(ii-1)*m:ii*m,1);
end

%
figure;
subplot(1,2,1)
hold on;
s = surf(bcX,bcY,uTrue);
s.EdgeColor = 'none';
surf(bcX,bcY,bcEdge,'EdgeColor','r');
title('True solution for u_{0,exact}')
xlabel('x-axis');
ylabel('y-axis');
set(gca,'fontsize',20);
subplot(1,2,2)
hold on;
s = surf(X,Y,UU);
s.EdgeColor = 'none';
surf(bcX,bcY,bcEdge,'EdgeColor','r');
title('Solution generated with 9-point Laplacian') 
xlabel('x-axis');
ylabel('y-axis');
set(gca,'fontsize',20);

%% Convergence for the 9-point Laplacian, fourth order Accuracy

a = 0; b = 1; c = 0; d = 1;
u0 = @(x,y) sin(4*pi*(x.^2 + y.^2)) + cos(4*pi*x.*y); % Define exact solution
fun = @(x,y) 16*pi*(cos(4*pi*(x.^2+y.^2))-pi*(x.^2+y.^2).*(4.*sin(4*pi*(x.^2+y.^2))+cos(4*pi*x.*y)));

for kk = 3:10
    
    h = 1/2^(kk); % grid point size
    hx(kk-2,1) = h;
    m = -1+1/h; % number of grid points
    n = m+2; % number of grid points with boundaries
    
    [X,Y] = meshgrid(a+h:h:b-h,c+h:h:d-h); % Define grid without boundaries
    [bcX,bcY] = meshgrid(a:h:b,c:h:d); % Define grid with boundaries 
    border = [ones(1,n);[ones(n-2,1),zeros(n-2,n-2),ones(n-2,1)];ones(1,n)]; % Define border
    uTrue = u0(bcX,bcY); % Define true solution

    fOld = zeros(m^2,1);
    for ii = 1:m
        fOld(1+(ii-1)*m:ii*m,1) = fun(X(ii,:),Y(ii,1));% + h^2/12*gun(X(ii,:),Y(ii,1));
    end
    LapRhs = LaplacianOfRhs(m,h,fun,bcX,bcY);
    f = fOld + h^2/12*LapRhs;%+h^2/12*fNew;

    bc = u0(bcX,bcY).*border;% + h^2/12* gun(bcX,bcY).*border; % Define baoundary at border
    
    A9 = poisson9(m);
    F = form_rhs9(m,f,bc);

    U = A9\F; % Calculate the solution U (as a vector)
    UU = zeros(m,m);
    for ii = 1:m % Reshape U as a m x m matrix
        UU(ii,:) = U(1+(ii-1)*m:ii*m,1);
    end
    
    GTE(kk-2,1) = max(max(abs(UU-uTrue(2:m+1,2:m+1)))); % Infinity norm 
end
disp('done')

%% Plot rate of convergence for 9-point Laplacian
figure;
subplot(1,2,1);
hold on;
plot(hx1,hx1.^2,'k-','linewidth',2)
plot(hx1,GTE1,'r--*','linewidth',2)
title('Rate of convergence of 9-point Laplacian')
ylabel('Global truncation error')
xlabel('stepsize, h')
legend('Expected Second order','Actual convergence','location','north')
grid on;
set(gca,'fontsize',20, 'XScale', 'log', 'YScale', 'log')
subplot(1,2,2);
hold on;
plot(hx,hx.^4,'k-','linewidth',2)
plot(hx,hx.^2,'b-','linewidth',2)
plot(hx,10^(-2).*GTE,'r--*','linewidth',2)
title('Rate of convergence for 9-point Laplacian, modified right hand side')
ylabel('Global truncation error')
xlabel('stepsize, h')
legend('Expected fourth order','Second order','Actual convergence','location','northwest')
grid on;
set(gca,'fontsize',20, 'XScale', 'log', 'YScale', 'log')
xlim([min(hx)*0.9 max(hx)*1.1])
