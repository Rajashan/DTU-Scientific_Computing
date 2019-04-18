function F = form_rhs9(m,f,bc) % Forming right hand side of AU=F
    h = 1/(m+1); % grid point size
    bc = bc/(6*h^2); % Normalize to prefactor of scheme
    F = f; % initialize F
    % Corner element at (x,y)=(0,0)
    F(1) = F(1) - bc(1,1)-4*bc(1,2)-bc(1,3)-4*bc(2,1)-bc(3,1); 
    % Corner element at (x,y)=(1,0)
    F(m) = F(m) - bc(1,m+2)-4*bc(1,m+1)-bc(1,m)-4*bc(2,m+2)-bc(3,m+2); 
    % Corner element at (x,y)=(0,1)
    F(1+m*(m-1)) = F(1+m*(m-1)) -...
                   bc(m+2,1)-4*bc(m+2,2)-bc(m+2,3)-4*bc(m+1,1)-bc(m,1); 
    % Corner element at (x,y)=(1,1)
    F(m^2) = F(m^2) -...
             bc(m+2,m+2)-4*bc(m+2,m+1)-bc(m+2,m)-4*bc(m+1,m+2)-bc(m,m+2); 
    % Edge elements at x=[0;1] and y=0
    F(2:m-1,1) = F(2:m-1,1) -... 
                 bc(1,2:m-1)'-4*bc(1,3:m)'-bc(1,4:m+1)'; 
    % Edge elements at x=0 and y=[0;1]
    F([m+1:m:1+(m-2)*m],1) = F([m+1:m:1+(m-2)*m],1) -... 
                             bc(2:m-1,1)-4*bc(3:m,1)-bc(4:m+1,1); 
    % Edge elements at x=1 and y=[0;1]
    F([2*m:m:m*(m-1)],1) = F([2*m:m:m*(m-1)],1) -... 
                           bc(2:m-1,m+2)-4*bc(3:m,m+2)-bc(4:m+1,m+2); 
    % Edge elements at x=[0;1] and y=1
    F(2+m*(m-1):(m-1)*(1+m)) = F(2+m*(m-1):(m-1)*(1+m)) -... 
                              bc(m+2,2:m-1)'-4*bc(m+2,3:m)'-bc(m+2,4:m+1)'; 
end