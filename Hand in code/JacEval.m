function J = JacEval(U,n,h,alpha,beta,epsilon) % Evaluate the Jacobian
    % First row (uses left boundary condition)
    J(1,1) = -2*epsilon/h^2 + (U(2)-alpha)/(2*h)-1;
    J(1,2) = epsilon/h^2 + U(1)/(2*h);
    for ii = 2:n-1 % All interior point
        J(ii,ii-1) = epsilon/h^2 - U(ii)/(2*h);
        J(ii,ii) = -2*epsilon/h^2 + (U(ii+1)-U(ii-1))/(2*h)-1;
        J(ii,ii+1) = epsilon/h^2 + U(ii)/(2*h);
    end
    % Last row (uses the right boundary condition)
    J(n,n-1) = epsilon/h^2 - U(n)/(2*h);
    J(n,n) = -2*epsilon/h^2 + (beta-U(n-1))/(2*h)-1;
end