function [x,Slow,divergent,k] = InexactNewtonsMethodODE(func, Jac, F, ii, tk,...
                                phi, dt, xinit, tol, minit, maxit, AT, varargin)
    % initialize
    k = 0;
    Slow = 0;
    divergent = 0;
    gamma = AT(2,2);
    x = xinit;
    % Compute the iteration matrix and LU factorize it
    J = feval(Jac,tk+dt,x,varargin{:}); % Evaluate the Jacobian
    M = eye(size(J))-dt*gamma.*J;
    [L,U] = lu(M,'vector');
    % Initial residual function
    F(:,ii) = feval(func,tk+dt,x,varargin{:});
    R = x - F(:,1:ii)*AT(1:ii,ii)*dt-phi;
    Convergent = 0;
    % Loop until convergent and below minimum and maximum number of 
    % Newton iterations
    while (~Convergent ||  k <= minit) && k < maxit
        k=k+1; % Coun the number of Newton iterations
        ROld = R; % Save the previous residual function
        dx = U\(L\R); % Calculate the change to the residual
        x = x - dx; % Change the x value
        F(:,ii) = feval(func,tk+dt,x,varargin{:}); % Evaluate function at new x
        R = x - F(:,1:ii)*AT(1:ii,ii)*dt - phi; % Calculate the residual
        alpha = norm(R,'inf')/norm(ROld,'inf'); % Estimate convergence rate
        if norm(R,'inf') < tol*0.1 % Test for convergence
            Convergent = 1;
        elseif alpha > 0.1 % Test for slow convergence
            disp('Slow !!!!')
            Slow = 1;
            break
        end
    end
    if k >= maxit % Test for divergence
        disp('Divergent !!!!')
        divergent = 1;
    end
end