function lim = dtlz_range(fname, M)
%DTLZ_RANGE Returns the decision range of a DTLZ function
%   The range is simply [0,1] for all variables. What varies is the number 
%   of decision variables in each problem. The equation for that is
%     n = (M-1) + k
%   wherein k = 5 for DTLZ1, 10 for DTLZ2-6, and 20 for DTLZ7.
%   
%   Syntax:
%      lim = get_range(fname, M)
%
%   Input arguments:
%      fname: a string with the name of the function ('dtlz1', 'dtlz2' etc.)
%      M: a scalar with the number of objectives
%
%   Output argument:
%      lim: a n x 2 matrix wherein the first column is the inferior limit 
%           (0), and the second column, the superior limit of search (1)

% Checks if the string has or not the prefix 'dtlz', or if the number later
% is greater than 7
fname = lower(fname);
if length(fname) < 5 || ~strcmp(fname(1:4), 'dtlz') || ...
   str2double(fname(5)) > 7
   error(['Sorry, the function ' fname ' is not implemented.'])
end

% If the name is o.k., defines the value of k
if strcmp(fname, 'dtlz1')
   k = 5;
elseif strcmp(fname, 'dtlz7')
   k = 20;
else %any other function
   k = 10;
end

n = (M-1) + k; %number of decision variables

lim = [zeros(n,1) ones(n,1)];

