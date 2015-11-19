# DEMO "Toolbox"
## Differential Evolution for Multiobjective Optimization

These codes were developed by [Fillipe Goulart](http://orcslab.cpdee.ufmg.br/index.php/current-members/46-fillipe-goulart-silva-mendes) ([fillipe.gsm@gmail.com](fillipe.gsm@gmail.com)) during his M.Sc. at Universidade Federal de Minas Gerais, under the mentoring of Prof. [Felipe Campelo](http://orcslab.cpdee.ufmg.br/index.php/faculty/5-felipe-campelo) ([fcampelo@ufmg.br](fcampelo@ufmg.br)). 

The _Octave-Matlab_ folder contains the implementations for Octave (which should work on Matlab too). The following algorithms are implemented:

- A posteriori methods (without preferences):  
    – DEMO [1]: the regular DEMO with non-dominated sorting;  
    – IBEA [2]: DEMO using indicators instead.  
- A priori or interactive (with preferences):  
    – R-DEMO [3]: R-NSGA-II but using the DEMO instead;  
    – PBEA [4]: IBEA but using a reference point;  
    – PAR-DEMO(nds) [5]: the method proposed by us, using nondominated sorting;  
    – PAR-DEMO(ε) [5]: the same method, but using indicators instead.  

Fillipe's M.Sc. thesis is available [here](http://ppgee.ufmg.br/defesas/1120M.PDF), and contains an extensive review on multiobjective optimization and preference-based methods. It also contains a more extensive description and discussion of the Preference-based Adaptive Region-of-interest (PAR) framework. 

If you use these codes in any way, please cite our paper [5]:

    @article{Goulart2016,
      doi = {10.1016/j.ins.2015.09.015},
      url = {http://dx.doi.org/10.1016/j.ins.2015.09.015},
      year  = {2016},
      month = {feb},
      publisher = {Elsevier {BV}},
      volume = {329},
      pages = {236--255},
      author = {Fillipe Goulart and Felipe Campelo},
      title = {Preference-guided evolutionary algorithms for many-objective optimization},
      journal = {Information Sciences}
    }



### References
1. T Robic and B Filipic. DEMO: Differential evolution for multiobjective optimization. Evolutionary Multi-Criterion Optimization, 520–533, 2005.  
1. Eckart Zitzler and S Kunzli. Indicator-based selection in multiobjective search. Parallel Problem Solving from Nature-PPSN VIII, (i):1–11, 2004.  
1. Kalyanmoy Deb, J. Sundar, Rao N. Udaya Bhaskara, and Shamik Chaudhuri. Reference Point Based Multi-Objective Optimization Using Evolutionary Algorithms. International Journal of Computational Intelligence Research, 2(3):273– 286, 2006.
1. Lothar Thiele, Kaisa Miettinen, PJ Korhonen, and Julian Molina. A preference- based evolutionary algorithm for multi-objective optimization. Evolutionary Computation, 17(3):411–436, 2009.
1. Fillipe Goulart and Felipe Campelo. Preference-guided evolutionary algorithms for many-objective optimization. Information Sciences, 329:236 – 255, 2016. Special issue on Discovery Science.
