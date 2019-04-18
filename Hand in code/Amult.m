function AU = Amult(U,m)
%Amult that employs natural ordering 

%Initial Values
h = 1/(m+1);
L = length(U);

%Construct Grid
I = zeros(m+2);

%Construct Inner Matrix
%UM = reshape(U,m,m)';
UM = rot90(reshape(U,m,m));

%Embed Inner Matrix
I(2:m+1,2:m+1) = UM;

%Construct Laplacian Approx (Negative Sign)
AU_Temp = (1/(h^2))*(4*UM - I(2:m+1,3:m+2) - I(2:m+1,1:m) - I(3:m+2,2:m+1) - I(1:m,2:m+1));


%Transform back to a vector
AU = reshape(rot90(AU_Temp,3),L,1);
end