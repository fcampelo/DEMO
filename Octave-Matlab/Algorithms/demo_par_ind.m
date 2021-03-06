function [fopt, xopt] = demo_par_ind(f, zr, xrange, options)
%DEMO_PAR_IND: Multi-objective optimization using the DEMO
%	Proposed method with Indicators.
%
%   Syntax:
%      [fopt, xopt] = demo_par_ind(f, zr, xrange)
%      [fopt, xopt] = demo_par_ind(f, zr, xrange, options)
%
%   Input arguments:
%      f: the objective function given as a handle, a .m file, inline, or
%         anything that can be computed using the "feval" function and can
%         handle multiple inputs. The output must be, for each point, a 
%         column vector of size m x 1, with m > 1 the number of objectives.
%      zr: a column vector with the aspiration point
%      xrange: a matrix with the inferior and superior limits of search.
%              If n is the dimension, it will be a n x 2 matrix, such
%              that the first column contains the inferior limit, and the
%              second, the superior one;
%      options: a struct with internal parameters of the algorithm:
%         .F: the scale factor to be used in the mutation (default: 0.5);
%         .CR: the crossover factor used in the recombination (def.: 0.3);
%         .mu: the population size (number of individuals) (def.: 100);
%         .kmax: maximum number of iterations (def.: 300);
%         .display: 'on' to display the population while the algorithm is
%                   being executed, and 'off' to not (default: 'off');
%         If any of the parameters is not set, the default ones are used
%         instead.
%
%   Output arguments:
%      fopt: the m x mu_opt matrix with the mu_opt best objectives
%      xopt: the n x mu_opt matrix with the mu_opt best individuals

% Check the parameters of the algorithm
if nargin < 4 %options was not provided
   options = struct();
end
options = check_input(options);

% Initial considerations
n = size(xrange,1); %dimension of the problem
P.x = rand(n, options.mu); %initial decision variables
P.f = fobjeval(f, P.x, xrange); %evaluates the initial population
m = size(P.f, 1); %number of objectives
k = 1; %iterations counter

% Beginning of the main loop
Pfirst = P;
while k <= options.kmax
   % Plot the current population (if desired)
   if strcmp(options.display, 'on')
      if m == 2
         plot(P.f(1,:), P.f(2,:), 'o', zr(1), zr(2), 'r*');
         title('Objective values during the execution')
         xlabel('f_1'), ylabel('f_2')
         drawnow
      elseif m == 3
         plot3(P.f(1,:), P.f(2,:), P.f(3,:), 'o', ...
         	zr(1), zr(2), zr(3), 'r*');
         title('Objective values during the execution')
         xlabel('f_1'), ylabel('f_2'), zlabel('f_3')
         drawnow
      end
   end
      
   % Perform the variation operation (mutation and recombination)
   O.x = mutation(P.x, options); %mutation
   O.x = recombination(P.x, O.x, options); %recombination
   O.x = repair(O.x); %assure the offspring do not cross the search limits
   O.f = fobjeval(f, O.x, xrange); %compute objective functions
   
   % Selection and updates
   P = selection(P, O, zr, options);
      
   fprintf('Iteration %d\n', k)
   k = k + 1;
end

% Return the final population
% First, unnormalize it
Xmin = repmat(xrange(:,1), 1, options.mu); %replicate inferior limit
Xmax = repmat(xrange(:,2), 1, options.mu); %replicate superior limit
Xun = (Xmax - Xmin).*P.x + Xmin;

% Then, return the nondominated front
ispar = ndset(P.f);
fopt = P.f(:,ispar);
xopt = Xun(:,ispar);

%=========================== Sub-functions ================================%
function phi = fobjeval(f, x, xrange)
%FOBJEVAL Evaluates the objective function
%   Since the population is normalized, this function unnormalizes it and
%   computes the objective values
%
%   Syntax:
%      phi = fobjeval(f, x, options)
%
%   Input arguments:
%      f: the objective function given as a handle, a .m file, inline, or
%         anything that can be computed using the "feval" function and can
%         handle multiple inputs
%      x: a n x mu matrix with mu individuals (points) and n variables 
%         (dimension size)
%      options: the struct with the parameters of the algorithm
%
%   Output argument:
%      phi: a m x mu matrix with the m objective values of the mu
%           individuals
 
mu = size(x, 2); %number of points
% Unnormalizes the population
Xmin = repmat(xrange(:,1), 1, mu); %replicates inferior limit
Xmax = repmat(xrange(:,2), 1, mu); %replicates superior limit
Xun = (Xmax - Xmin).*x + Xmin;

phi = feval(f, Xun);
%--------------------------------------------------------------------------%
function Xo = mutation(Xp, options)
%MUTATION Performs mutation in the individuals
%   The mutation is one of the operators responsible for random changes in
%   the individuals. Each parent x will have a new individual, called trial
%   vector u, after the mutation.
%   To do that, pick up two random individuals from the population, x2 and
%   x3, and creates a difference vector v = x2 - x3. Then, chooses another
%   point, called base vector, xb, and creates the trial vector by
%
%      u = xb + F*v = xb + F*(x2 - x3)
%
%   wherein F is an internal parameter, called scale factor.
%
%   Syntax:
%      Xo = mutation(Xp, options)
%
%   Input arguments:
%      Xp: a n x mu matrix with mu "parents" and of dimension n
%      options: the struct with the internal parameters
%
%   Output arguments:
%      Xo: a n x mu matrix with the mu mutated individuals (of dimension n)

% Creates a mu x mu matrix of 1:n elements on each row
A = repmat((1:options.mu), options.mu, 1);
% Now, as taken at the MatLab Central, one removes the diagonal of
% A, because it contains indexes that repeat the current i-th
% individual
A = A';
A(logical(eye(size(A)))) = []; %removes the diagonal
A = transpose(reshape(A, options.mu-1, options.mu)); %reshapes

% Now, creates a matrix that permutes the elements of A randomly
[~, J] = sort(rand(size(A)),2);
Ilin = bsxfun(@plus,(J-1)*options.mu,(1:options.mu)');
A(:) = A(Ilin);

% Chooses three random points (for each row)
xbase = Xp(:, A(:,1)); %base vectors
v = Xp(:,A(:,2)) - Xp(:,A(:,3)); %difference vector

% Performs the mutation
Xo = xbase + options.F*v;
%--------------------------------------------------------------------------%
function Xo = recombination(Xp, Xm, options)
%RECOMBINATION Performs recombination in the individuals
%   The recombination combines the information of the parents and the
%   mutated individuals (also called "trial vectors") to create the
%   offspring. Assuming x represents the i-th parent, and u the i-th trial
%   vector (obtained from the mutation), the offspring xo will have the
%   following j-th coordinate:
%
%      xo_j = u_j if rand_j <= CR
%             x_j otherwise
%
%   wherein rand_j is a number drawn from a uniform distribution from 0 to
%   1, and CR is called the crossover factor. To prevent mere copies, at
%   least one coordinate is guaranteed to belong to the trial vector.
%
%   Syntax:
%      Xo = recombination(Xp, Xm, options)
%
%   Input arguments:
%      Xp: a n x mu matrix with the mu parents
%      Xm: a n x mu matrix with the mu mutated points
%      options: the struct with the internal parameters
%
%   Output argument:
%      Xo: a n x mu matrix with the recombinated points (offspring)

% Draws random numbers and checks whether they are smaller or
% greater than CR
n = size(Xp, 1); %dimension of the problem
aux = rand(n, options.mu) <= options.CR;
% Now assures at least one coordinate will be changed, that is,
% there is at least one 'true' in each column
auxs = sum(aux) == 0; %gets the columns with no trues
indc = find(auxs); %get the number of the columns
indr = randi(n, 1, sum(auxs)); %define random indexes of rows
if isempty(indr), indr = []; end
if isempty(indc), indc = []; end
ind = sub2ind([n, options.mu], indr, indc); %converts to indexes
aux(ind) = true;

% Finally, creates the offspring
Xo = Xp;
Xo(aux) = Xm(aux);
%--------------------------------------------------------------------------%
function Xo = repair(Xo)
%REPAIR Truncates the population to be in the feasible region
%
%   Syntax:
%      Xo = repair(Xo, options)

% This is easy, because the population must be inside the interval [0, 1]
Xo = max(Xo, 0); %corrects inferior limit
Xo = min(Xo, 1); %superior limit

%--------------------------------------------------------------------------%
function Pnew = selection(P, O, zr, options)
%SELECTION Selects the next population
%   Each parent is compared to its offspring. If the parent dominates its 
%   child, then it goes to the next population. If the offspring dominates 
%   the parent, that new member is added. However, if they are incomparable
%   (there is no mutual domination), them both are sent to the next 
%   population. After that, the new set of individuals must be truncated to 
%   mu, wherein mu is the original number of points.
%   This is accomplished by the use of "non-dominated sorting", that is,
%   ranks the individual in fronts of non-domination, and within each
%   front, measures them by using crowding distance. With regard to these
%   two metrics, the best individuals are kept in the new population.
%
%   Syntax:
%      Pnew = selection(P, O, options)
%
%   Input arguments:
%      P: a struct with the parents (x and f)
%      O: a struct with the offspring
%      options: the struct with the algorithm's parameters
%
%   Output argument:
%      Pnew: the new population (a struct with x and f)

% ------ First part: checks dominance between parents and offspring
% Verifies whether parent dominates offspring
aux1 = all(P.f <= O.f, 1);
aux2 = any(P.f < O.f, 1);
auxp = and(aux1, aux2); %P dominates O
% Now, where offspring dominates parent
aux1 = all(P.f >= O.f, 1);
aux2 = any(P.f > O.f, 1);
auxo = and(aux1, aux2); %O dominates P
auxpo = and(~auxp, ~auxo); %P and O are incomparable
% New population (where P dominates O, O dominates P and where they are 
% incomparable)
R.f = [P.f(:,auxp), O.f(:,auxo), P.f(:,auxpo), O.f(:,auxpo)];
R.x = [P.x(:,auxp), O.x(:,auxo), P.x(:,auxpo), O.x(:,auxpo)];

% ------- Second part: get the points belonging to the ROI
[Fn, zrn] = fstandardize(R.f, zr, options.pareto_ranges); %standardize population
[iroi, s] = find_roi(Fn, zrn, options.rho); %get the ROI
Rroi.f = R.f(:,iroi);
Rroi.x = R.x(:,iroi);

% ------- Third part: corrects the population size
% If the number of points inside the ROI is mu, stop.
% If there are too few points, pick more according to their ASF value;
% Finally, if there are too many, correct the size with the Non-dominated 
% sorting
mu_roi = size(Rroi.f, 2);
if mu_roi == options.mu %right number of points (this will happen so rarely...)
	Pnew = Rroi; %both .f and .x

elseif mu_roi < options.mu %too few points
	% Grab the points outside the ROI
	R.f = R.f(:,~iroi); R.x = R.x(:,~iroi); s = s(~iroi);
	[~, is] = sort(s);
	muaux = options.mu - mu_roi; %missing number of points
	Pnew.f = [Rroi.f, R.f(:,is(1:muaux))];
	Pnew.x = [Rroi.x, R.x(:,is(1:muaux))];

else %too many points
	% Truncate the population using the indicator procedure
	musofar = size(Rroi.f,2); %number of points so far
	for ii = 1:(musofar - options.mu)
	   F_ind = eps_fitness(Rroi.f, options.kappa);
	   [~, imin] = min(F_ind);
	   % Removes this worst point
	   Rroi.f(:,imin) = [];
	   Rroi.x(:,imin) = [];
	end
	
	Pnew = Rroi; %both .f and .x
end

%--------------------------------------------------------------------------%
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

% Finally, computes the value F by summing the corresponding columns
F = sum(-exp(-I_eps/kappa), 2);

%--------------------------------------------------------------------------%
function options = check_input(options)
%CHECK_INPUT Checks the parameters of the algorithm before
%   This sub-function checks the endogenous parameters of the algorithm. If
%   they are not set, the default ones are used

if ~isfield(options, 'F') %scale factor
   options.F = 0.5;
end

if ~isfield(options, 'CR') %crossover factor
   options.CR = 0.3;
end

if ~isfield(options, 'kmax') %maximum number of iterations
   options.kmax = 300;
end

if ~isfield(options, 'mu') %population size
   options.mu = 100;
end

if ~isfield(options, 'display') %show or not the population during execution
   options.display = 'off';
end

if ~isfield(options, 'rho') %augmenting factor for the ASF
	options.rho = 1e-6;
end

if ~isfield(options, 'kappa') %fitness scaling factor for the indicator
   options.kappa = 0.05; %a random number from my head
end

% If the ranges of the Pareto front are given (ideal and Nadir points), use them
% Otherwise, the default will be an empty matrix, and the algorithm will use the
% current population to approximate them.
if ~isfield(options, 'pareto_ranges')
	options.pareto_ranges = [];
end
