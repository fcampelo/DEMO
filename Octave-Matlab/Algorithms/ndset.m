function ispar = ndset(F)
%NDSET Finds the nondominated set of a set of objective points
%   
%   Syntax:
%      ispar = ndset(F)
%
%   Input argument:
%      F: a m x mu matrix with mu points and m objectives
%
%   Output argument:
%      ispar: a mu-length vector with true in the nondominated points

mu = size(F,2); %number of points

% The idea is to compare each point with the other ones
f1 = permute(F, [1 3 2]); %puts in the 3D direction
f1 = f1(:,ones(mu,1),:);
f2 = F(:,:,ones(mu,1));

% Now, for the ii-th slice, the ii-th individual is compared with all of the
% others at once. Then, the usual operations of domination are checked
% Checks where f1 dominates f2
aux1 = all(f1 <= f2,1);
aux2 = any(f1 < f2,1);
auxf1 = and(aux1, aux2);
% Checks where f1 is dominated by f2
aux1 = all(f1 >= f2,1);
aux2 = any(f1 > f2,1);
auxf2 = and(aux1, aux2);

% dom will be a 3D matrix (1 x mu x mu) such that, for the ii-th slice, it
% will contain +1 if fii dominates the current point, -1 if it is dominated 
% by it, and 0 if they are incomparable
dom = zeros([1 mu mu]);
dom(auxf1) = 1;
dom(auxf2) = -1;

% Finally, the slices with no -1 are nondominated
ispar = all(dom ~= -1, 2);
ispar = ispar(:);
