function [iroi, s] = find_roi(phi, zr, rho)
%FIND_ROI Finds the points that belong to the ROI
%    Given a set of objective values phi and a reference point zr, this
%    function finds the points that belong to the region of interest (ROI)

mu = size(phi,2);
if nargin < 3, rho = 1e6; end

s = asf(phi, zr, rho); %computes the ASF
[~, imin] = min(s); %finds the best value
d = dist2zr(phi, zr);
dmin = d(imin);
% Computes the new reference points by adding dmin to each coordinate
% (remember we are in the objective space)
m = size(zr,1); %number of objectives
zraux = repmat(zr, 1, m);
zraux = zraux + diag(dmin(ones(m,1)))*sign(s(imin));
% Now, finds the closest points to each new reference point. The point with 
% highest s will give the maximum admissible asf to form the ROI
d = asf(phi, zraux, rho);
[~, imin] = min(d, [], 2);
fl = phi(:,imin); %fl contains the limiting points

fmin = min(fl, [], 2); fmax = max(fl, [], 2);
Fmin = fmin(:, ones(mu,1)); Fmax = fmax(:, ones(mu,1)); %replications
aux1 = any(phi < Fmin, 1);
aux2 = any(phi > Fmax, 1);
iroi = ~or(aux1, aux2);


%==========================================================================%
function d = dist2zr(F, zr)
%DIST2ZR Computes the distance to each point to the reference point(s)
%   This was taken from Matlab Tips and Tricks

% One needs to transpose the matrix so it works with the code
F = F'; zr = zr';
m = size(F,1); n = size(zr,1);
F = permute(F, [1 3 2]);
zr = permute(zr, [3 1 2]);
d = sqrt(sum(abs(F(:, ones(1, n), :) - zr(ones(1, m), :, :)).^2, 3));
d = d';
