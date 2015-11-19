function f = dtlz2(x, M)
%DTZL2 DTLZ2 multi-objective function
%   This function represents a hyper-sphere.
%   Using k = 10, the number of dimensions must be n = (M - 1) + k.
%   The Pareto optimal solutions are obtained when the last k variables of x
%   are equal to 0.5.
%
%   Syntax:
%      f = dtzl2(x, M)
%
%   Input arguments:
%      x: a n x mu matrix with mu points and n dimensions
%      M: a scalar with the number of objectives
%
%   Output argument:
%      f: a m x mu matrix with mu points and their m objectives computed at
%         the input

k = 10;
% Error check: the number of dimensions must be M-1+k
n = (M-1) + k; %this is the default
if size(x,1) ~= n
   error(['Using k = 10, it is required that the number of dimensions be'...
   ' n = (M - 1) + k = %d in this case.'], n)
end

xm = x(n-k+1:end,:); %xm contains the last k variables
g = sum((xm - 0.5).^2, 1);

% Computes the functions
f(1,:) = (1 + g).*prod(cos(pi/2*x(1:M-1,:)),1);
for ii = 2:M-1
   f(ii,:) = (1 + g) .* prod(cos(pi/2*x(1:M-ii,:)),1) .* ...
      sin(pi/2*x(M-ii+1,:));
end
f(M,:) = (1 + g).*sin(pi/2*x(1,:));
