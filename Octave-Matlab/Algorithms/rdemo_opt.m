function [fopt, xopt] = rdemo_opt(f, zr, xrange, options)
%RDEMO_OPT: Multi-objective optimization using the DEB method
%
%   Syntax:
%      [fopt, xopt] = rdemo_opt(f, zr, xrange)
%      [fopt, xopt] = rdemo_opt(f, zr, xrange, options)
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
      
#   warning('Iteration %d', k)
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

#warning('==== The Global Search finished with %d iterations ===.', k)

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

% ------- Second part: get the distance of the points to zr
% When I say distance, it is really a distance, an Euclidean one!
[Fn, zrn] = fstandardize(R.f, zr, options.pareto_ranges); %standardize population
aux = bsxfun(@minus, Fn, zrn);
d = sqrt( sum( aux.^2, 1) );

% ------- Third part: truncates the population according to distance
% I am not in the mood to make this part more beautiful... So, this works, and 
% since I do not intent to use this algorithm for anything more important, I 
% will just keep it like this.
Daux = d; %makes no sense, but I need it like this...
Faux = R.f;
Xaux = R.x;
keepit = true;
Fnew = [];
Dnew = [];
Xnew = [];
musofar = 0; %individuals in the new population so far

while keepit
   ispar = ndset(Faux); %gets the current 
   if sum(ispar) + musofar < options.mu %the population is not complete
      Fnew = [Fnew, Faux(:,ispar)];
      Dnew = [Dnew, Daux(:,ispar)];
      Xnew = [Xnew, Xaux(:,ispar)];
      % Removes the current nondominated set from the population
      Faux(:,ispar) = [];
      Daux(:,ispar) = [];
      Xaux(:,ispar) = [];
      musofar = size(Fnew,2);      
   else %there are too many individuals in the current front
      % From this front, gets points according to the epsilon-dominance, and
      % complements the population with them
      Faux = Faux(:,ispar); Daux = Daux(:,ispar); Xaux = Xaux(:,ispar);
      while true                  
         iseps = eps_dominance(Faux, Daux, options);            
         if sum(iseps) < options.mu - musofar %this is not enough
            Fnew = [Fnew, Faux(:,iseps)]; 
            Dnew = [Dnew, Daux(:,iseps)];
            Xnew = [Xnew, Xaux(:,iseps)];
            Faux(:,iseps) = []; Daux(:,iseps) = []; Xaux(:,iseps) = [];
            musofar = size(Fnew,2);            
         else % there is more than enough points to complement the new P
            % In that case, chooses the ones with the least distance
            Faux = Faux(:,iseps); Daux = Daux(:,iseps); Xaux = Xaux(:,iseps);
            [~, isort] = sort(Daux);
            muaux = options.mu - musofar;
            Fnew = [Fnew, Faux(:,isort(1:muaux))];
            Dnew = [Dnew, Daux(:,isort(1:muaux))];          
            Xnew = [Xnew, Xaux(:,isort(1:muaux))];
            break;
         end
      end      
      keepit = false;
   end
end 
     
% Updates the population
Pnew.f = Fnew;
Pnew.x = Xnew;

% ---------------------------------------------------------------------------- %
function iseps = eps_dominance(Fnd, Dnd, options)
%EPS_DOMINANCE Applies the "epsilon dominance" in the current front
%   Remember it is assumed that the incoming individuals are all 
%   nondominated.

% Finds the boxes of each point, and compares the ii-th box with all of the
% others. This will be done at once using 3D matrices
mu = size(Fnd, 2); %number of points in this front
b = floor(Fnd./options.epsilon)*options.epsilon;
b1 = permute(b, [1 3 2]); %rotates to the 3D direction
b1 = b1(:, ones(mu,1), :); %replicates columnwise
b2 = b(:,:,ones(mu,1)); %replicates in the 3D direction

% Now, checks the individuals within the same box: auxbox will be a matrix
% such that the ii-th row contains true in the individuals that are in same
% box as the individual ii (notice that will be at least a diagonal matrix)
auxbox = all(b1 == b2,1);
auxbox = permute(auxbox, [3 2 1]);

d = Dnd(ones(mu,1),:); %replicates the distances for mu individuals
d(~auxbox) = inf;

% From each box, gets the points with the least distance
[~, imin] = min(d, [], 2);

% Finally, imin contains, for each row, the number of the point that
% has the minimum distance to the reference point. If, for example, in the
% row 2 imin contains the number 3, then it means that the individual 3 is
% inside the same box, and it has the smallest distance. So, the point 2 is
% not optimal in this case
aux = (1:mu)';
iseps = (imin == aux);

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


if ~isfield(options, 'epsilon') %variable to define the DM's ROI
   options.epsilon = 1e-3;
end

% If the ranges of the Pareto front are given (ideal and Nadir points), use them
% Otherwise, the default will be an empty matrix, and the algorithm will use the
% current population to approximate them.
if ~isfield(options, 'pareto_ranges')
	options.pareto_ranges = [];
end
