function G = GEval(u,m,h,alpha,beta,epsilon)
    % Apply equation (2.106) in the book
    % First point (uses the left boundary condition)
    G(1) = epsilon*(alpha-2*u(1)+u(2))/h^2 + u(1)*((u(2)-alpha)/(2*h)-1);
    for ii = 2:m-1 % All interior points
        G(ii) = epsilon*(u(ii-1)-2*u(ii)+u(ii+1))/h^2 + u(ii)*((u(ii+1)-u(ii-1))/(2*h)-1);
    end
    % Last point (uses the right boundary condition)
    G(m) = epsilon*(u(m-1)-2*u(m)+beta)/h^2 + u(m)*((beta-u(m-1))/(2*h)-1);
    G=G';
end
