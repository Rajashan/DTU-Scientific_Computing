function LapRhs = LaplacianOfRhs(m,h,fun,bcX,bcY) % 5-point Laplacian of right hand function
    Lap5fun = zeros(m,m); % Initialize (m x m)-matrix
    for ii = 1:m % loop over rows (from bottom to top)
    for jj = 1:m % loop over columbs (from right to left)
        % Evaluate fun (right hand function) using the 5-point Lapacian scheme
        Lap5fun(ii,jj) = (fun(bcX(ii,jj+1),bcY(ii,jj+1))+fun(bcX(ii+1,jj),bcY(ii+1,jj))+...
                         fun(bcX(ii+1,jj+2),bcY(ii+1,jj+2))+fun(bcX(ii+2,jj+1),bcY(ii+2,jj+1))-...
                         4*fun(bcX(ii+1,jj+1),bcY(ii+1,jj+1)))/h^2;
    end
    end
    LapRhs = zeros(m^2,1); % Initialize (m^2 x 1)-vector as result
    for ii = 1:m % Reshape to a (m^1 x 1)-vector using the natural rowwise order
        LapRhs(1+(ii-1)*m:ii*m,1) = Lap5fun(ii,:);
    end
end