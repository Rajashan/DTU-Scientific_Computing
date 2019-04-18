function A = poisson9(m)
    e = ones(m,1);
    S = spdiags([-e -10*e -e], [-1 0 1], m, m);
    I = spdiags([-1/2*e e -1/2*e], [-1 0 1], m, m);
    A = 1/6*(m+1)^2*(kron(I,S)+kron(S,I));
end