function g = VanderpolDiffusion1(t,x,p)
% State independent diffusion for the Van der pol equation. 

sigma = p(2);
g = [0.0;sigma];