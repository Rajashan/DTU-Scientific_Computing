function [tout,yout,stats] = ESDIRK23solver(func,Jac,ta,tb,h0,x0,...
                                  abstol,reltol,eps,facmin,facmax,varargin)
    % ESDIRK23 solver with adaptive step size
    s = 3; p = 2;
    % Define the Butcher tableau
    gamma = (2-sqrt(2))/2;
    A = [[0,0,0];...
         [gamma,gamma,0];...
         [sqrt(2)/4,sqrt(2)/4,gamma]];
    AT = transpose(A);
    b = [sqrt(2)/4;sqrt(2)/4;gamma];
    bhat = [(4-sqrt(2))/12;(4+3*sqrt(2))/12;(2-sqrt(2))/6];
    c = [0;2*gamma;1];
    d = b-bhat;
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
    % use the ESDIRK23 method
    m = 0;
    while t < tb
        m=m+1; % Count number of steps
        stats.fEval(m,1) = 0;
        stats.JEval(m,1) = 0;
        if t+h>tb
            h = tb-t;
        end
        n = 0;
        Accept = false;
        while ~Accept % Run until error is below tolerance
            n = n+1; % Count number of rejected steps (-1)
            Slow = 0; divergent = 0;
            % Stage 1
            T(1) = t; 
            X(:,1) = x;
            F(:,1) = feval(func,T(1),X(:,1),varargin{:});
            % Stage 2 and 3
            T(2:s) = t + h*c(2:s);
            for ii = 2:s
                F(:,ii) = F(:,ii-1);
                xinit = X(:,ii-1) + F(:,1:ii)*AT(1:ii,ii).*h; % Initial guess
                % Call inexact Newton solver
                [X(:,ii),Slow,divergent,k] = InexactNewtonsMethodODE(func, Jac, F,...
                    ii , t, x, h*c(ii), xinit, abstol, 5, 1000, AT, varargin{:});
                F(:,ii) = feval(func,T(ii),X(:,ii),varargin{:}); % Evaluated the function
                stats.k(m,ii-1) = k; % Save the number of Newton interations
                if Slow || divergent % Check for divergence or slow convergence
                    break
                end
            end
            stats.fEval(m,1) = stats.fEval(m,1) + 1 + 2*(1+1+k);
            stats.JEval(m,1) = stats.JEval(m,1) + 2;
            if ~Slow && ~divergent % If not slow or divergent
                x = X(:,3); % save the next value for the solution
                error = F*d*h; % estimate the error
                r = max(abs(error)./max(abstol,abs(x).*reltol)); % Calculate the r value
                Accept = (r <= 1); % Test we r<=1, within tolerance
            end
            if Accept % If within tolerance
                % To print progress
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
                % Progress to the next time step
                t = t+h;
                % Svae stat values
                stats.n(end+1,1) = n;
                stats.R = [stats.R;r];
                stats.tau = [stats.tau;error'];
                stats.hx = [stats.hx;h];
                % Save to the output
                tout = [tout;t];
                yout = [yout;x'];
                % Chanve the step size
                if n == 1
                    h = max([facmin,min([facmax,(eps/r).^(1/(p+1))])])*h;
                else
                    h = max([facmin,min([facmax,...
                        (h/stats.hx(end-1)).*(eps*stats.R(end-1)/r^2).^(1/(p+1))])])*h;
                end
            elseif Slow || divergent
                % Decrease step size of slow congerce or divergence in
                % Newton solver
                h = h*facmin;
            else
                % lower step size if error is to large
                h = max([facmin,min([facmax,(eps/r).^(1/(p+1))])])*h;
            end
        end
    end
end