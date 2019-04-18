clear all;
close all;
%% 2-point BVPs: Newton's method for solving nonlinear BVPs
alpha = -1; % Left boundary condition
beta = 1.5; % Right boundary condition
epsilon = 0.1; % Value for epsilon factor (epsilon << 1)
a = 0; % Starting point of interval
b = 1; % End point of interval
m = 100; % Number of grid points
h=(b-a)/(m+1); % Step size
X = a+h:h:b-h; % Grid
w0 = 1/2*(a-b+beta-alpha); % Width of the internal layer (omega_0)
xbar = 1/2*(a+b-beta-alpha); % Position of the center of the internal layer
maxit = 1000; % Maximum number of newton iterations
tol = 10^(-3); % Tolerance of Newton solver

uOuter(:,1) = [a,X,b]' + alpha - a; % Solution with only left boundary condition
uOuter(:,2) = [a,X,b]' + beta - b; % Solution with only right boundary condition
f = @(x) x - xbar + w0*tanh(w0*(x-xbar)/(2*epsilon)); % Initial guess (eq. 2.105)
U = f(X)'; % Set initial guess

% Newton solver
G = GEval(U,m,h,alpha,beta,epsilon) % Generate the residual for the system
J = JacEval(U,m,h,alpha,beta,epsilon); % Define the Jacobian
k=0; % Counter of Newton iterations
while ((norm(G,'inf') > tol) && (k <= maxit)) % Newton solver
    k=k+1; % Count Newton step
    dGdx = J; % Estimate the change in gradient of G
    dx = -dGdx\G; % Estimate the change in G
    U = U + dx; % Change the solution U
    J = JacEval(U,m,h,alpha,beta,epsilon); % Estimate the Jacobian
    G = GEval(U,m,h,alpha,beta,epsilon); % reevaluate the residual
    if k == maxit % Test if reached maximum number of iterations
        disp('Used max iterations !!!!!!')
    end
end
U=[alpha;U;beta]; % Add the boundary conditions to the solution

%% Plot the results
figure;
hold on;
plot([a,X,b],uOuter(:,1),'k--','linewidth',2);
plot([a,X,b],uOuter(:,2),'k-.','linewidth',2);
plot([a X b],U,'r-*','linewidth',3)
fplot(f,[a b],'b-','linewidth',2)
title(strcat('Solution with \epsilon=',num2str(epsilon)))
legend('Solution at left boundary','Solution at right boundary',...
    strcat('Approximated solution (k=',num2str(k),')'),'Initial guess','location','northwest')
xlabel('x')
ylim([min(U) max(f(X))*1.2])
set(gca,'fontsize',20)



