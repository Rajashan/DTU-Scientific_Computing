clear all; close all; clc

% Exact Solution and RHS
u = @(x,y) exp(pi*x).*sin(pi*y)+0.5*(x.*y).^2;   %True Solution
f = @(x,y) x.^2 + y.^2;                          %Right Hand Side
m = 2^6 - 1;                                     %Grid
U = zeros(m*m,1);
F = -form_rhs(f,u,m);

% Define region
a = 0; b = 1;
c = 0; d = 1;

%Grid Points and Grid
n = m+2;        % Number of grid points with boundaries
h = 1/(m+1);    % Grid point size
[X,Y] = meshgrid(a+h:h:b-h,c+h:h:d-h); % Define interior grid
[bcX,bcY] = meshgrid(a:h:b,c:h:d);     % Define boundaries' grid

%Solution and Values in the Grid
u = @(x,y) exp(pi*x).*sin(pi*y)+0.5*(x.*y).^2;           %True Solution
%uTrue = u0(bcX,bcY);                                    % Define true solution
border = [ones(1,n);[ones(n-2,1),zeros(n-2,n-2),ones(n-2,1)];ones(1,n)]; % Define border
bc = u(bcX,bcY).*border;                               % Define baoundary at border
bcEdge = bc;                                            % Define the edge as the boundary
bcEdge(bc == 0) = NaN;                                  % Remove every other point than the boundary

%Form right hand side
fun = @(x,y) x.^2 + y.^2;                %Right Hand Side
f = zeros(m^2,1);                        % Initialize f as a m^2 x 1 vector
for ii = 1:m                             % Create the value for the f-vector
    f(1+(ii-1)*m:ii*m,1) = fun(X(ii,:),Y(ii,1));
end
%(Negative sign of RHS)
%F = -form_rhs5(m,f,bc);

%{
%Plot Graph
h=1/(m+1);
[X,Y] = meshgrid(0:h:1,0:h:1);
figure(1)
surf(X,Y,u(X,Y))
%}

omega = 0.5;
%omega = 1/3;
epsilon = 1.0E-10;

for i = 1:5
    R = F + Amult(U,m);
    fprintf('*** Outer iteration: %3d, rel. resid.: %e\n', ...
        i, norm(R,2)/norm(F,2));
    if (norm(R,2)/norm(F,2) < epsilon)
        break;
    end
    U = Vcycle2(U,omega,3,m,F);
    %plotU(m,U,u);
    pause(.5);
end



function plotU(m,U)
h=1/(m+1);
x=linspace(1/h,1-1/h,m);
y=linspace(1/h,1-1/h,m);
[X,Y]=meshgrid(x,y);
surf(X, Y, reshape(U,[m,m])');
shading interp;
title('Computed solution');
xlabel('x');
ylabel('y');
zlabel('U');
end