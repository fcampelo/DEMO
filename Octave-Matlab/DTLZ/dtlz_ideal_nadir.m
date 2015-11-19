function pareto_ranges = dtlz_ideal_nadir(fname, M)
%DTLZ_IDEAL_NADIR Get the true limits of the Pareto front
%	 This function returns the ideal and the Nadir points of a DTLZ function.
%	 For the DTLZ1:
%	 	 zstar = [0 0 0 ... 0]
%		 znad = [0.5 0.5 0.5 ... 0.5]
%	
%	 For the DTLZ2 through DTLZ4:
%		 zstar = [0 0 0 ... 0]
%		 znad = [1 1 1 ... 1]
%
%	 DTLZ5 and DTLZ6: the Nadir is a little bit complicated
%		 zstar = [0 0 0 ... 0]
%		 znad = (1/sqrt(2)) ^ [M-2, M-2, M-3, ..., 3, 2, 1, 0]
%
%	 DTLZ7: This is the most trickier.
%	 The ideal is composed of M-1 zeros, and the last component is equal to 
%		 faux = 2*(M - (M-1)*xaux/2*(1 + sin(3*pi*xaux))), 
%		 	 wherein xaux = 0.85940, obtained numerically by me
%	 The Nadir has M-1 components equal to xaux (not faux), and the last one is 
%	 2*M. So,
%		 zstar = [0 0 0 ... 0 faux]
%		 znad = [xaux xaux xaux ... xaux 2*M]
%
%	 Syntax:
%		 pareto_ranges = dtlz_ideal_nadir(fname, M)
%


switch (fname)
	
	case 'dtlz1'
		zstar = zeros(M,1);
		znad = 0.5(ones(M,1));
	
	case {'dtlz2', 'dtlz3', 'dtlz4'}
		zstar = zeros(M,1);
		znad = ones(M,1);
		
	case {'dtlz5', 'dtlz6'}
		zstar = zeros(M,1);
		znad = (1/sqrt(2)) .^ ([M-2, M-2:-1:0]');
	
	case 'dtlz7'
		xaux = 0.85940;
		faux = 2*(M - (M-1)*xaux/2*(1 + sin(3*pi*xaux)));
		zstar = [zeros(M-1,1); faux];
		znad = [xaux(ones(M-1,1)); 2*M];
	
	otherwise
		error('Sorry, what the heck of a function is %s?', fname)
end

pareto_ranges = [zstar, znad];
