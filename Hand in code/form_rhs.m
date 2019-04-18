function Fi = form_rhs(fun,ufun,m)

%Construct Full Grid
h = 1/(m+1);
[Xc,Yc] = meshgrid(0:h:1,0:h:1);

%Compute U (Everywhere but only boundary needed)
F = ufun(Xc,Yc);

%Compute F in interior points and transform into column vector
Fi = fun(Xc(2:m+1,2:m+1),Yc(2:m+1,2:m+1));
Fi = reshape(rot90(Fi,3),m^2,1);

%Correct All Corners
Fi(1) = Fi(1) - (1/h^2)*(F(m+1,1) + F(m+2,2));      %LowerLeft Corner = First Element
Fi(m) = Fi(m) - (1/h^2)*(F(m+2,m+1) + F(m+1,m+2));  %LowerRight Corner = M'th Element
M_ur = m^2-m+1;
Fi(M_ur) = Fi(M_ur) - (1/h^2)*(F(1,2) + F(2,1));    %UpperLeft Corner = (M^2 - M + 1)'th element
Fi(m^2) = Fi(m^2) - (1/h^2)*(F(1,m+1) + F(2,m+2));  %UpperRight Corner = M^2'th element

%Correct All Non-Corners
Fi(2:m-1) = Fi(2:m-1) - (1/h^2)*F(m+2,3:m)';                        %Bottom Corrected
Fi(M_ur+1:m^2-1) = Fi(M_ur+1:m^2-1) - (1/h^2)*F(1,3:m)';            %Upper Corrected
Fi(1+m:m:M_ur-m) = Fi(1+m:m:M_ur-m) - (1/h^2)*fliplr((F(3:m,1)));   %Left Corrected
Fi((2*m:m:m^2-m)) = Fi(2*m:m:m^2-m) - (1/h^2)*fliplr((F(3:m,m+2))); %Right Corrected
end

%Fi(1) = Fi(1) - (1/h^2)*(F(2,1) + F(1,2)); % steffens guess
%Fi(m) = Fi(m) - (1/h^2)*(F(1,m+1) + F(2,m+2)); % Steffens guess
%Fi(M_ur) = Fi(M_ur) - (1/h^2)*(F(m+2,2) + F(m+1,1)); % Steffens guess
%Fi(m^2) = Fi(m^2) - (1/h^2)*(F(m+2,m+1) + F(m+1,m+2)); % Steffens guess
%Fi(1:m) = Fi(1:m) - (1/h^2)*F(1,2:m+1)';          %Bottom, steffen
%Fi(M_ur:m^2) = Fi(M_ur:m^2) - (1/h^2)*(F(m+2,2:m+1))'; %Upper, steffen
% Fi(1+m:m:(m^2-2*m+1)) = Fi(1+m:m:(m^2-2*m+1)) - (1/h^2)*fliplr((F(3:m,1))); % Left, steffen