%% Solve -A*U = -F with PCG (Preconditioned Conjugate Gradient Method)

%Initialize. Function Handles and Interior Points
m = 3;
fun = @(x,y) 16*pi.*(cos(4*pi*(x.^2+y.^2))-pi.*(x.^2+y.^2).*(4.*sin(4*pi*(x.^2+y.^2))+cos(4*pi*x.*y)));
u0 = @(x,y) sin(4*pi*(x.^2 + y.^2)) + cos(4*pi*x.*y);

%Construct the Right-Hand Side
F = -form_rhs(fun,u0,m);

%Solve the Problem using PCG
U = zeros(1,9);                 %Initial Guess
J = @(U) Amult(U,m);            %Function Handle
[U_Sol,~,~,Iter,~] = pcg(J,F);  %Call PCG-Matlab Function

