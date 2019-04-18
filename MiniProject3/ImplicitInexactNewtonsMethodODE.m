function [x,Slow,divergent,k] = ImplicitInexactNewtonsMethodODE(func, Jac,...
                        T, dt, xinit, tol, minit, maxit, s, A, cc, varargin)
% initialize
    k = 0;
    Slow = 0;
    divergent = 0; 
    % Shape x and c vectors
    x = []; c= [];
    for ii = 1:s
        x = [x;xinit];
        c = [c;cc(ii);cc(ii)];
    end
    phi = x;
    % Compute the iteration matrix and LU factorize it
    J = feval(Jac,T(s),x,varargin{:});
    M = zeros(size(J)*s);
    AI = zeros(size(J)*s);
    % Shape the M and A matrix
    for ii = 1:s
        M(2*ii-1:2*ii,2*ii-1:2*ii) = eye(size(J))-dt.*J;
        for jj = 1:s
            AI(2*ii-1:2*ii,2*jj-1:2*jj) = A(ii,jj).*eye(size(J));
        end
    end
    AIT = transpose(AI);
    % LU factorize
    [L,U,pp] = lu(M,'vector');
    for ii = 1:s % shape F vector
        F(1,2*ii-1:2*ii) = feval(func,T(ii),xinit,varargin{:});
    end
    x = x + F*c*dt; % Backward Eulor for first guess
    for ii = 1:s % shape F vector
        F(1,2*ii-1:2*ii) = feval(func,T(ii),x(2*ii-1:2*ii),varargin{:});
    end
    % Initialize the residual function
    R = x - transpose(F*AIT).*dt-phi;
    Convergent = 0;
    % Loop until converged or/and minit < k < maxit
    while (~Convergent ||  k <= minit) && k < maxit
        k=k+1;
        ROld = R; % Save preveous residual
        dx = U\(L\R); % Estimate change in x
        x = x - dx; % Change x
        for ii = 1:s % shape F vector
            F(1,2*ii-1:2*ii) = feval(func,T(ii),x(2*ii-1:2*ii),varargin{:});
        end
        R = x - transpose(F*AIT).*dt - phi; % update residual
        alpha = norm(R,'inf')/norm(ROld,'inf'); % Measure rate of convergence
        if norm(R,'inf') < tol*0.1 % Test id converged
            Convergent = 1;
        elseif alpha > 0.1 % Test if slow convergence
            Slow = 1;
            break
        end
    end
    if k >= maxit % Test if diverged
        disp('Divergent !!!!!')
        divergent = 1;
    end
end