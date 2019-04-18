function F = form_rhs5(m,f,bc)
    h = 1/(m+1); % grid point size
    bc = bc/h^2; % Rescale the boundary conditions
    F = f; % Initialize right hand side without boundaries
    F(1) = F(1) - bc(1,2) - bc(2,1); % add boundary in left lower corner
    F(m) = F(m) - bc(m+1,1) - bc(m+2,2); % add boundary in right lower corner
    for ii = 2:m-1 % add the boundaries at the bottom
        F(ii) = F(ii) - bc(1,ii+1);
    end
    for ii = 2:m-1 % add the boundaries at the left and right sides
        F(1+(ii-1)*m) = F(1+(ii-1)*m) - bc(ii+1,1); % left side
        F(ii*m) = F(ii*m) - bc(ii+1,m+2); % right side
    end
    F(1+m*(m-1)) = F(1+m*(m-1)) - bc(m+1,1) - bc(m+2,2); % add boundary at left upper corner
    F(m^2) = F(m^2) - bc(m+2,m+1) - bc(m+1,m+2); % add boundary at right upper corner
    %for ii = 2:m % add the boundaries at the upper side
    for ii = 2:(m-1) %Correction Phillip
        F(ii+m*(m-1)) = F(ii+m*(m-1)) - bc(m+2,ii+1);
        %Der er fejl her, der trækkes F(5,4) fra to gange. Dvs F(m^2) er for lille // Phillip
    end
end