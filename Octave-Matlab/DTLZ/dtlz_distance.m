function d = dtlz_distance(xopt, fname)
%DTLZ_DISTANCE Finds the distance of a point to the optimum
%   Since, in the DTLZ functions, the last k variables of a decision vector
%   represent the "distance" to the Pareto set and the first ones the
%   "position" in the front, a good way to find the distance to the optimum
%   is to check these last variables.
%   For each function, the last k variables must be a given value, like 0.5,
%   1, 0 etc. Hence, the distance from xopt(end-k+1:end) to 0.5, 0, whatever
%   will give how far the points are to the optimal set.
%
%   Syntax:
%      d = dtlz_distance(xopt, fname)
%
%   Input arguments:
%      xopt: a n x mu matrix of mu individuals of dimension n (it must be
%            in accordance with the required dimension given by Deb)
%      fname: a string with the name of the function (e.g., 'dtlz1')
%
%   Output argument:
%      d: a mu-vector with the distance of each point to the optimal set

% For each function, a different value of k and reference point is given.
% So, I will just use a stupid switch and everybody's happy
mu = size(xopt,2); %number of points
switch fname
   
   case 'dtlz1' %k = 5 and x(last) = 0.5
      k = 5;
      xlast = repmat(0.5, [k mu]);
      d = sum((xopt(end-k+1:end,:) - xlast).^2, 1);
      
   case 'dtlz2' %k = 10 and x(last) = 0.5
      k = 10;
      xlast = repmat(0.5, [k mu]);
      d = sum((xopt(end-k+1:end,:) - xlast).^2, 1);
      
   case 'dtlz3' %k = 10 and x(last) = 0.5
      k = 10;
      xlast = repmat(0.5, [k mu]);
      d = sum((xopt(end-k+1:end,:) - xlast).^2, 1);
   
   case 'dtlz4' %k = 10 and x(last) = 0.5
      k = 10;
      xlast = repmat(0.5, [k mu]);
      d = sum((xopt(end-k+1:end,:) - xlast).^2, 1);
      
   case 'dtlz5' %k = 10 and x(last) = 0.5
      k = 10;
      xlast = repmat(0.5, [k mu]);
      d = sum((xopt(end-k+1:end,:) - xlast).^2, 1);
      
   case 'dtlz6' %k = 10 and x(last) = 0
      k = 10;
      xlast = zeros([k mu]);
      d = sum((xopt(end-k+1:end,:) - xlast).^2, 1);
      
   case 'dtlz7' %k = 20 and x(last) = 0
      k = 20;
      xlast = zeros([k mu]);
      d = sum((xopt(end-k+1:end,:) - xlast).^2, 1);
end
