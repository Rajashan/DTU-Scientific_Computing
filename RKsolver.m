function [tout,yout] = RKsolver( func, tspan, n, y0)
    t = linspace(tspan(1),tspan(2),n); % time interval
    h = t(2)-t(1); % find stepsize
    % initialize outputs
    tout = zeros(n,1);
    yout = zeros(n,1);
    % First value is set to initial value
    yout(1) = y0;
    tout(1) = t(1);
    % use the Runge Kutta method to get up to n values
    for ii = 2:n
        s1 = func(t(ii-1),yout(ii-1));
        s2 = func(t(ii-1)+1/2*h,yout(ii-1)+1/2*h*s1);
        s3 = func(t(ii-1)+1/2*h,yout(ii-1)+1/2*h*s2);
        s4 = func(t(ii-1)+h,yout(ii-1)+h*s3);
        yout(ii) = yout(ii-1)+h/6*(s1+2*s2+2*s3+s4); % set y-value
        tout(ii) = t(ii); % set time-value
    end
end