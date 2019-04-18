% problem size and spacing
m = 100;
h = 1/(m+1);
lambda = zeros(m,m);


% Nested for-loops evaluating the eigenvalues for the 2D-Poisson problem
for p = 1:1:m
    for q = 1:1:m
        lambda(p,q) = 1+(h^2/4)*(2/h^2).*((cos(p*pi*h)-1)+(cos(q*pi*h)-1));
    end
end

% Only considering the higher frequency eigenvalues
lambda2 = lambda(m/2+1:end,m/2+1:end);

% We take omega in [0,2] and preallocate
omega = 0:0.01:2;
index = zeros(1,length(omega));

% Maximum absolute eigenvalue for each value of omega. 
i = 1;
for k = omega
    temp=max(abs((1-k).*eye(length(lambda2))+k.*lambda2));
    index(1,i)=max(temp);
    i = i+1;
end

figure(1)
plot(omega,index)
xlabel('omega')
ylabel('max(gamma_{pq})')