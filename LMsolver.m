function [tout,yout] = LMsolver( func, tspan, n, y0)
    t = linspace(tspan(1),tspan(2),n);% time interval
    h = t(2)-t(1); % find setpsize

    alpha = [-1,0,1]; % define alphs coefficients
    beta = h.*[1/3, 4/3, 1/3]; % defince beta coefficients
    % initialize outputs
    tout = zeros(n,1);
    yout = zeros(n,1);
    % First value is set to initial value
    yout(1) = y0;
    tout(1) = t(1);
    % use the Runge Kutta method to get second value
    s1 = func(t(1),yout(1));
    s2 = func(t(1)+1/2*h,yout(1)+1/2*h*s1);
    s3 = func(t(1)+1/2*h,yout(1)+1/2*h*s2);
    s4 = func(t(1)+h,yout(1)+h*s3);
    yout(2) = yout(1)+h/6*(s1+2*s2+2*s3+s4);
    tout(2) = t(2);
    % Use method from part d) to get up to n values
    for ii = 3:n
        known = alpha(1)*yout(ii-2)-beta(2)*func(t(ii-1),yout(ii-1))-beta(1)*func(t(ii-2),yout(ii-2));
        tosolve = @(z) alpha(3)*z -beta(3)*func(t(ii),z) + known;
        z0 = [0 -2*known];
        o = -2*known;
        try % try to find a solution to the implicit system
            z = fzero(tosolve,z0); % use fzero to find solution to step value
            yout(ii) = z;
            tout(ii) = t(ii);
        catch % if none out off solution and break
            yout = yout(1:ii-1);
            tout = tout(1:ii-1);
            break
        end 
    end
end