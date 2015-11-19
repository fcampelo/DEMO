function cdist = crowdingdistance(F)
%CROWDINGDISTANCE Computes the crowding distance of a nondominated front
%   The crowding distance gives a measure of how close the individuals are
%   with regard to its neighbors. The higher this value, the greater the
%   spacing. This is used to promote better diversity in the population.
%
%   Prototype:
%      cdist = crowdingdistance(F)
%
%   Input arguments:
%      F: an m x mu matrix with mu individuals and m objectives
%
%   Output argument:
%      cdist: a m-length column vector

[m, mu] = size(F); %gets the size of F

if mu == 2
   cdist = [inf inf]';
   return
end

[Fs, Is] = sort(F,2); %sorts the objectives by individuals
% Creates the numerator
C = Fs(:,3:end) - Fs(:,1:end-2);
C = [inf(m,1) C inf(m,1)]; %complements with inf in the extremes

% Indexing to permute the C matrix in the right ordering
Aux = repmat((1:m)', 1, mu);
ind = sub2ind([m, mu], Aux(:), Is(:));
C2(ind) = C(:);
C = reshape(C2, [m, mu]);

% Constructs the denominator
den = repmat(Fs(:,end) - Fs(:,1), 1, mu);

% Calculates the crowding distance
cdist = sum(C./den,1);
cdist = cdist(:); %assures a column vector
