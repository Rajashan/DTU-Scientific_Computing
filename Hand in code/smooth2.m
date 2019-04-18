function Unew = smooth2(U,omega,m,F)

h = 1/(m+1);
Unew = U - omega * ( (h^2/4)*Amult(U,m) + (h^2/4)*F);

end