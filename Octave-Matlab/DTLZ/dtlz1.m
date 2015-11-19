function f = dtlz1(x, M)
%DTZL1 DTLZ1 multi-objective function
%   This function represents a hyper-plane.
%   Using k = 5, the number of dimensions must be n = (M - 1) + k.
%   The Pareto front of this function happens at xm := x(n-k+1:end) = 0.5, 
%   that is, when the last k variables are all equal to 0.5. The first 
%   M-1 variables are varied to map the whole Pareto front. Of course, this 
%   is not the easiest thing to do for more than M = 3 objectives. However, 
%   the hyper-surface is a hyper-plane satisfying the equation
%      f1 + f2 + f3 + ... + fm = 0.5
%
%   Syntax:
%      f = dtzl1(x, M)
%
%   Input arguments:
%      x: a n x mu matrix with mu points and n dimensions
%      M: a scalar with the number of objectives
%
%   Output argument:
%      f: a m x mu matrix with mu points and their m objectives computed at
%         the input

k = 5; %as suggested by Deb
% Error check: the number of dimensions must be M-1+k
n = (M-1) + k; %this is the default
if size(x,1) ~= n
   error(['Using k = 5, it is required that the number of dimensions be '...
   'n = (M - 1) + k = %d in this case.'], n)
end

xm = x(n-k+1:end,:); %xm contains the last k variables
g = 100*(k + sum((xm - 0.5).^2 - cos(20*pi*(xm - 0.5)),1));

% Now, computes the functions. The first and the last will be
% written separately to facilitate things
f(1,:) = 1/2*prod(x(1:M-1,:),1).*(1 + g);
for ii = 2:M-1
   f(ii,:) = 1/2*prod(x(1:M-ii,:),1).*(1 - x(M-ii+1,:)).*(1 + g);
end
f(M,:) = 1/2*(1 - x(1,:)).*(1 + g);
