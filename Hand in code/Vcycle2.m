function Unew = Vcycle2(U,omega,nsmooth,m,F)
% Approximately solve: A*U = F

l2m = log2(m+1);
assert(l2m==round(l2m));    %m must be m = 2^k - 1 otherwise error
assert(length(U)==m*m);     %u must have length m^2
if (m==1)
    % if we are at the coarsest level
    % TODO: solve the only remaining equation directly!
        Unew = F/16;
else
    % 1. TODO: pre-smooth the error
    %    perform <nsmooth> Jacobi iterations
        for n = 1:nsmooth
            U = smooth2(U,omega,m,F);
        end
    % 2. TODO: calculate the residual
        R = - (Amult(U,m) + F);
    % 3. TODO: coarsen the residual
        Rc = coarsen(R,m);
    % 4. recurse to Vcycle on a coarser grid
        mc = (m-1)/2;
        Ecoarse = Vcycle2(zeros(mc*mc,1),omega,nsmooth,mc,-Rc);
    % 5. TODO: interpolate the error
        %E = Interpolate(Ecoarse,mc);
        E = interpolate_a(Ecoarse,2*mc+1);
    % 6. TODO: update the solution given the interpolated error
        U = U + E;
    % 7. TODO: post-smooth the error
    %    perform <nsmooth> Jacobi iterations
        for n = 1:nsmooth
            Unew = smooth2(U,omega,m,F);
        end
end
end