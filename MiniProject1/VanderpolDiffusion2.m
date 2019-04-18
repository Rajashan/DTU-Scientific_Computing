function g = VanderpolDiffusion2(t,x,p)
% State dependent diffusion for the Van der pol equation. 

sigma = p(2);
g = [0.0; sigma*(1.0+x(1)*x(1))];