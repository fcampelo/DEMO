function f = dtlz7(x, M)
%DTZL7 DTLZ7 multi-objective function
%   This function has a disconnected Pareto front
%   Using k = 20, the number of dimensions must be n = (M - 1) + k.
%   The Pareto optimal solutions are obtained when the last k variables of x
%   are equal to 0.
%
%   Syntax:
%      f = dtzl7(x, M)
%
%   Input arguments:
%      x: a n x mu matrix with mu points and n dimensions
%      M: a scalar with the number of objectives
%
%   Output argument:
%      f: a m x mu matrix with mu points and their m objectives computed at
%         the input

k = 20;
% Error check: the number of dimensions must be M-1+k
n = (M-1) + k; %this is the default
if size(x,1) ~= n
   error(['Using k = 20, it is required that the number of dimensions be'...
   ' n = (M - 1) + k = %d in this case.'], n)
end

% Writes down the auxiliar function g
xm = x(n-k+1:end,:); %xm contains the last k variables
g = 1 + 9/k*sum(xm,1);

% Now, computes the first M-1 objective functions
f(1:M-1,:) = x(1:M-1,:);

% The last function requires another auxiliar variable
gaux = g(ones(M-1,1),:); %replicates the g function
h = M - sum(f./(1+gaux).*(1 + sin(3*pi*f)),1);
f(M,:) = (1 + g).*h;
