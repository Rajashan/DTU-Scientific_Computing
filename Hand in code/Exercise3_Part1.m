%Exercise 3 Part1

% Define region
a = 0; b = 1;
c = 0; d = 1;

%Grid Points and Grid
m = 3;          % Number of interior grid points
n = m+2;        % Number of grid points with boundaries
h = 1/(m+1);    % Grid point size
[X,Y] = meshgrid(a+h:h:b-h,c+h:h:d-h); % Define interior grid
[bcX,bcY] = meshgrid(a:h:b,c:h:d);     % Define boundaries' grid

%Solution and Values in the Grid
u0 = @(x,y) sin(4*pi*(x.^2 + y.^2)) + cos(4*pi*x.*y);    % Define exact solution
border = [ones(1,n);[ones(n-2,1),zeros(n-2,n-2),ones(n-2,1)];ones(1,n)]; % Define border
bc = u0(bcX,bcY).*border;                               % Define baoundary at border
bcEdge = bc;                                            % Define the edge as the boundary
bcEdge(bc == 0) = NaN;                                  % Remove every other point than the boundary

%Form right hand side
fun = @(x,y) 16*pi.*(cos(4*pi*(x.^2+y.^2))-pi.*(x.^2+y.^2).*(4.*sin(4*pi*(x.^2+y.^2))+cos(4*pi*x.*y)));
f = zeros(m^2,1);       % Initialize f as a m^2 x 1 vector
for ii = 1:m            % Create the value for the f-vector
    f(1+(ii-1)*m:ii*m,1) = fun(X(ii,:),Y(ii,1));
end
%(Negative sign of RHS)
F1 = -form_rhs5(m,f,bc);

% m = 3;
% fun = @(x,y) 16*pi.*(cos(4*pi*(x.^2+y.^2))-pi.*(x.^2+y.^2).*(4.*sin(4*pi*(x.^2+y.^2))+cos(4*pi*x.*y)));
F = -form_rhs(fun,u0,m);

%Solve -A*U = -F with PCG (Preconditioned Conjugate Gradient Method)
U = zeros(1,9);         %Initial Guess
J = @(U) Amult(U,m);    %Function Handle
U_sol1 = pcg(J,F);       %Call PCG-Matlab Function
J = @(U) Amult2(U,m);    %Function Handle
U_sol2 = pcg(J,F1);

[F1 F]
[U_sol1 U_sol2]

%%


%% Interpolation Test
m = 3;
h = 1/(m+1);

[X,Y] = meshgrid(0+h:h:1-h,0+h:h:1-h);
Z = sin(X*pi)-sin(Y*pi);

[R,mf] = Interpolate(Z,m);
z = rot90(reshape(R,mf,mf));

%{
figure;
subplot(1,2,1)
surf(X,Y,Z)
title('True solution')
subplot(1,2,2)
surf(Xf,Yf,z)
title('Approximated solution')
%}

%% Shit Works
