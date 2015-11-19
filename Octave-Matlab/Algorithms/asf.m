function s = asf(phi, zr, rho)
%ASF Achievement Scalarizing Function
%   This function receives the multi-objective function and converts it
%   into a scalar function by means of the ASF, given by
%
%      s(f(x), zr) = max[(f(x) - zr)] + rho*sum[(f(x) - zr)]
%
%   Remember: this only makes sense when the objectives (and the reference 
%   point) are standardized, which is assumed in this function.
%
%   Syntax:
%      s = asf(phi, zr)
%      s = asf(phi, zr, rho)
%
%   Input Arguments:
%	   phi: a m x mu matrix with the objectives of the population
%      zr: a m x nr matrix with nr reference points. They can be
%          feasible, infeasible or anything;
%      rho: the augmenting factor (def.: 1e-6)
%
%   Output arguments:
%      s: a 1 x mu vector with the ASF value for each of the mu individuals
%
%	See also FSTANDARDIZE.

[m, mu] = size(phi); %number of objectives and of points
nr = size(zr, 2); %number of reference points

% First, gets some default values
if nargin < 3, rho = 1e-6; end

% This function is prepared to the fact that more than one reference point
% is adopted. Then, one puts each reference point in a 3D slice, and
% replicates it to the number of points. Also, the objectives are repeated
% in the 3D direction for each reference point
zr = permute(zr, [1 3 2]);
zr = repmat(zr, [1 mu 1]);
phi = repmat(phi, [1 1 nr]);

% First, infinite norm
sinf = max(phi - zr, [], 1);

% Second, the linear term to assure properly Pareto optimal solutions
s1 = rho* sum(phi - zr, 1);

% Finally, combines everything and permutes to become a 2D matrix, wherein
% the ii-th column contains the ASF value for the ii-th individual in each
% reference point
s = sinf + s1;
s = permute(s, [3 2 1]);
