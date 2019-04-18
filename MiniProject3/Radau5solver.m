function [tout,yout,stats] = Radau5solver(func,Jac,ta,tb,h0,x0,...
                                  abstol,reltol,eps,facmin,facmax,varargin)
    % Radau5 solver with adaptive step size
    s = 3; p = 5; % Number of stages and order of method
    A = [[(88-7*sqrt(6))/360,(296-169*sqrt(6))/1800,(-2+3*sqrt(6))/225];...
        [(296+169*sqrt(6))/1800,(88+7*sqrt(6))/360,(-2-3*sqrt(6))/225];...
        [(16-sqrt(6))/36,(16+sqrt(6))/36,1/9]]; % Butcher Tableau A matrix
    c = [(4-sqrt(6))/10;(4+sqrt(6))/10;1]; % Butcher Tableau c matrix
    % Initial values
    t = ta;
    h = h0;
    x = x0;
    % Initial value stats output
    stats.n = [];
    stats.k = [];
    stats.R = 1;
    stats.tau = x'.*0;
    stats.hx = h0;
    stats.fEval = [];
    stats.JEval = [];
    % To print progress
    tstep = (tb-ta)*0.1; tn = 1;
    % First value is set to initial value
    yout = x';
    tout = t;
    m=0;
    while t < tb
        m=m+1; % Count number of steps
        stats.fEval(m,1) = 0;
        stats.JEval(m,1) = 0;
        if t+h>tb
            h = tb-t;
        end
        n = 0;
        Accept = false;
        while ~Accept
            n = n+1; % Number of rejected time steps
            % inner time steps
            T(1:s) = t + h*c(1:s);
            % Call the implicit inexact Newton solver
            [X,Slow,divergent,k] = ImplicitInexactNewtonsMethodODE(func,...
                Jac, T, h, x, min(abstol, reltol), 10, 1000, s, A, c, varargin{:});
            X = reshape(X,length(x),s); % Reshape output
            stats.k(m,1) = k; % Save number of Newton iterations
            x = X(:,3); % Set next solution value
            stats.fEval(m,1) = stats.fEval(m,1) + 6 + 3*k;
            stats.JEval(m,1) = stats.JEval(m,1) + 1;
            % If Newton method converged with a acceptable rate
            if (~divergent && ~Slow)
                % Use step doubling
                hm = 0.5*h; % Half step size
                Tm = t + hm.*c; % New inner time steps
                % Call the implicit inexact Newton solver
                [Xhat,~,~,kk] = ImplicitInexactNewtonsMethodODE(func,...
                    Jac, Tm, hm, x, min(abstol, reltol), 10, 1000, s, A, c, varargin{:});
                stats.fEval(m,1) = stats.fEval(m,1) + 6 + 3*kk;
                stats.JEval(m,1) = stats.JEval(m,1) + 1;
                Xhat = reshape(Xhat,length(x),s);
                xhat = Xhat(:,3); % Save solution from step doubling
                error = xhat - x; % Estimate error
                % Test if within toleraces
                r = max(abs(error)./max(abstol,abs(x).*reltol)); 
                Accept = (r <= 1); 
            end
            if Accept % If accepted step
                % Print progress
                if t >= tstep
                    if tn == 1
                        fprintf(strcat('\n',num2str(tn)))
                    elseif tn == 9
                        fprintf(strcat(num2str(tn),'\n'))
                    else
                        fprintf(num2str(tn))
                    end
                    tn = tn+1;
                    tstep = tstep + (tb-ta)*0.1;
                end
                t = t+h; % Set next time step
                % Save stats
                stats.n(end+1,1) = n;
                stats.R = [stats.R;r];
                stats.tau = [stats.tau;error'];
                stats.hx = [stats.hx;h];
                tout = [tout;t];
                yout = [yout;x'];
                % PI control change in step size
                if n == 1
                    h = max([facmin,min([facmax,(eps/r).^(1/(p+1))])])*h;
                else
                    h = max([facmin,min([facmax,...
                        (h/stats.hx(end-1)).*(eps*stats.R(end-1)/r^2).^(1/(p+1))])])*h;
                end
            elseif divergent || Slow % If diverged or slow convergence rate
                h = h*facmin; % Reduce step size
            else % If step did not get accepted
                h = max([facmin,min([facmax,(eps/r).^(1/(p+1))])])*h;
            end
        end
    end
end
    
    
    