function F = eps_fitness(phi, kappa)
%EPS_FITNESS Computes fitness according to I_epsilon indicator
%   The fitness for the point x_i is computed as
%
%      F(x_i) = sum_{x_j, j ~= i} -exp(-I(x_j, x_i)/kappa)
%
%   where in I(x_j, x_i) is the minimum value such that x_j dominates x_i
%   in the "weakly" sense; kappa is a fitness scaling such that no one 
%   knows exactly a good value for it.
%
%   This function performs an adaptative version, wherein the objectives are
%   normalized to the interval [0, 1], and so the indicators. In this way,
%   the algorithm becomes not so dependent on the parameter kappa.

mu = size(phi,2); %number of points

% Finds the upper and lower bounds for each objective, and replicates them
% with the number of points. Then, normalize the objectives
fmax = max(phi, [], 2); fmax = fmax(:,ones(mu,1)); %upper bound
fmin = min(phi, [], 2); fmin = fmin(:,ones(mu,1)); %lower bound
phi = (phi - fmin)./(fmax - fmin); %normalization

% Replicates phi to allow multiple comparisons at once
f1 = phi(:,:,ones(mu,1)); %replicates in the 3D direction
f2 = permute(phi, [1 3 2]); %rotates in to the 3D direction
f2 = f2(:,ones(mu,1),:); %replicates

% Computes the epsilon indicator, finds the greatest absolute value, and 
% normalizes I_eps by it
I_eps = max(f1 - f2, [], 1); %gets the required I_epsilon 
I_eps = permute(I_eps, [3 2 1]); %turns into a 2D matrix
c = max(abs(I_eps), [], 2); %finds the maximum for each individual
I_eps = I_eps./c(:,ones(mu,1));

% Now, removes the diagonals, because it is an indicator with respect to 
% itself, and it will always be zero
I_eps = I_eps';
I_eps(logical(eye([mu mu]))) = []; %removes the main diagonal
I_eps = transpose(reshape(I_eps, [mu-1 mu])); %reshapes

I_eps,

% Finally, computes the value F by summing the corresponding columns
F = sum(-exp(-I_eps/kappa), 2);