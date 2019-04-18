function [xmean,s,xmeanp2s,xmeanm2s]=ScalarSampleMeanStdVar(x)

% Calculation of mean, std and pm 2*std at each instance, hence dimension 3
xmean = mean(x,3);
s = std(x,0,3);
xmeanp2s = xmean + 2*s;
xmeanm2s = xmean - 2*s;