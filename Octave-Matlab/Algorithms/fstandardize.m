function [Fn, zrn] = fstandardize(phi, zr, pareto_ranges)
%FSTANDARDIZE Standardize the objectives
%	 Gets the ideal and Nadir points of a problem in order to standardize the 
%	 objectives (and the reference point). If the input pareto_ranges is not 
%	 provided or empty then use the following approximations:
%
%		 z*_app = min(phi, [], 2) - 0.05
%		 z^nad_app = max(phi, [], 2)
%
%	 The standardization is given by
%	 	 Fn_i = (phi_i - z*)/(znad - z*), i = 1, 2, ..., mu

if nargin < 3, pareto_ranges = []; end

% Check if the limits of the Pareto front were given
if isempty(pareto_ranges)
	ispar = ndset(phi); %get the non-dominated front
	zstar = min(phi(:,ispar), [], 2) - 0.05;
	znad = max(phi(:,ispar), [], 2);
else
	zstar = pareto_ranges(:,1);
	znad = pareto_ranges(:,2);
end

% Now, standardize the population and the reference point
zaux = znad - zstar;
zrn = (zr - zstar)./zaux;
Fn = bsxfun(@minus, phi, zstar)./repmat(zaux, 1, size(phi,2));
	
