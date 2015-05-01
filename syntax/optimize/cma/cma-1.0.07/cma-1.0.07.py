   1  #!/usr/bin/env python 
   2  """Module cma implements the CMA-ES (Covariance Matrix Adaptation Evolution Strategy). 
   3   
   4  CMA-ES is a stochastic optimizer for robust non-linear non-convex 
   5  derivative- and function-value-free numerical optimization. 
   6   
   7  This implementation can be used with Python versions 2.6, 2.7, 3.x. 
   8   
   9  CMA-ES searches for a minimizer (a solution x in :math:`R^n`) of an 
  10  objective function f (cost function), such that f(x) is minimal. 
  11  Regarding f, only a passably reliable ranking of the candidate 
  12  solutions in each iteration is necessary. Neither the function values 
  13  itself, nor the gradient of f need to be available or do matter (like 
  14  in the downhill simplex Nelder-Mead algorithm). Some termination 
  15  criteria however depend on actual f-values. 
  16   
  17  Two interfaces are provided: 
  18   
  19    - function `fmin(func, x0, sigma0,...)` 
  20          runs a complete minimization 
  21          of the objective function func with CMA-ES. 
  22   
  23    - class `CMAEvolutionStrategy` 
  24        allows for minimization such that the 
  25        control of the iteration loop remains with the user. 
  26   
  27   
  28  Used packages: 
  29   
  30      - unavoidable: `numpy` (see `barecmaes2.py` if `numpy` is not 
  31        available), 
  32      - avoidable with small changes: `time`, `sys` 
  33      - optional: `matplotlib.pyplot` (for `plot` etc., highly 
  34        recommended), `pprint` (pretty print), `pickle` (in class 
  35        `Sections`), `doctest`, `inspect`, `pygsl` (never by default) 
  36   
  37  Install 
  38  ------- 
  39  The file ``cma.py`` only needs to be visible in the python path (e.g. in 
  40  the current working directory), but can also be installed in the 
  41  terminal command line by:: 
  42   
  43      python cma.py --install 
  44   
  45  which solely calls the ``setup`` function from the standard 
  46  ``distutils.core`` package for installation. If the ``setup.py`` 
  47  file is been provided with ``cma.py``, the standard call is 
  48   
  49      python setup.py cma 
  50   
  51  Both calls need to see ``cma.py`` in the current working directory and 
  52  might need to be preceded with ``sudo``. 
  53   
  54  We can install or upgrade the currently installed version also with:: 
  55   
  56      pip install --upgrade cma 
  57   
  58  Testing 
  59  ------- 
  60  From the system shell:: 
  61   
  62      python cma.py --test 
  63   
  64  or from the Python shell ``ipython -pylab``:: 
  65   
  66      run cma.py --test 
  67   
  68  or from any python shell 
  69   
  70      import cma 
  71      cma.main('--test') 
  72   
  73  runs ``doctest.testmod(cma)`` showing only exceptions (and not the 
  74  tests that fail due to small differences in the output) and should 
  75  run without complaints in about between 20 and 100 seconds. 
  76   
  77  Example 
  78  ------- 
  79  From a python shell:: 
  80   
  81      import cma 
  82      help(cma)  # "this" help message, use cma? in ipython 
  83      help(cma.fmin) 
  84      help(cma.CMAEvolutionStrategy) 
  85      help(cma.CMAOptions) 
  86      cma.CMAOptions('tol')  # display 'tolerance' termination options 
  87      cma.CMAOptions('verb') # display verbosity options 
  88      res = cma.fmin(cma.Fcts.tablet, 15 * [1], 1) 
  89      res[0]  # best evaluated solution 
  90      res[5]  # mean solution, presumably better with noise 
  91   
  92  :See: `fmin()`, `CMAOptions`, `CMAEvolutionStrategy` 
  93   
  94  :Author: Nikolaus Hansen, 2008-2014 
  95   
  96  :License: MIT, see below. 
  97   
  98  """ 
  99   
 100  # The MIT License (MIT) 
 101  # Copyright (c) 2014 Inria 
 102  # Author: Nikolaus Hansen, 2008-2014 
 103  # 
 104  # Permission is hereby granted, free of charge, to any person obtaining 
 105  # a copy of this software and associated documentation files (the 
 106  # "Software"), to deal in the Software without restriction, including 
 107  # without limitation the rights to use, copy, modify, merge, publish, 
 108  # distribute, sublicense, and/or sell copies of the Software, and to 
 109  # permit persons to whom the Software is furnished to do so, subject to 
 110  # the following conditions: 
 111  # 
 112  # The above copyright and authorship notice and this permission notice 
 113  # shall be included in all copies or substantial portions of the 
 114  # Software. 
 115  # 
 116  # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
 117  # EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
 118  # MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT 
 119  # . IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR 
 120  # ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
 121  # CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION 
 122  # WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 123   
 124   
 125  # (note to self) for testing: 
 126  #   pyflakes cma.py   # finds bugs by static analysis 
 127  #   pychecker --limit 60 cma.py  # also executes, all 60 warnings checked 
 128  #   or python ~/Downloads/pychecker-0.8.19/pychecker/checker.py cma.py 
 129  #   python cma.py -t -quiet # executes implemented tests based on doctest 
 130  #   python -3 cma.py --test  2> out2to3warnings.txt # 
 131   
 132  # to create a html documentation file: 
 133  #    pydoc -w cma  # edit the header (remove local pointers) 
 134  #    epydoc cma.py  # comes close to javadoc but does not find the 
 135  #                   # links of function references etc 
 136  #    doxygen needs @package cma as first line in the module docstring 
 137  #       some things like class attributes are not interpreted correctly 
 138  #    sphinx: doc style of doc.python.org, could not make it work (yet) 
 139   
 140  # TODO: should optimize return self? Disadvantage: getting a hand on 
 141  #       the logger reference becomes very obscure in a one-line call. 
 142  # TODO: separate display and logging options, those CMAEvolutionStrategy 
 143  #       instances don't use them themselves (probably all) 
 144  # TODO: disp method is implemented in CMAEvolutionStrategy and in 
 145  #       CMADataLogger separately, OOOptimizer.disp_str should return a str 
 146  #       which can be used uniformly? 
 147  # TODO: check scitools.easyviz and how big the adaptation would be 
 148  # TODO: split tell into a variable transformation part and the "pure" 
 149  #       functionality 
 150  #       usecase: es.tell_geno(X, [func(es.pheno(x)) for x in X]) 
 151  #       genotypic repair is not part of tell_geno 
 152  # TODO: copy_always optional parameter does not make much sense, 
 153  #       as one can always copy the input argument first, 
 154  #       however some calls are simpler 
 155  # TODO: generalize input logger in optimize() as after_iteration_handler 
 156  #       (which is logger.add by default)? One difficulty is that 
 157  #       the logger object is returned (not anymore when return of optimize 
 158  #       is change). Another difficulty is the obscure usage of modulo 
 159  #       for writing a final data line in optimize. 
 160  # TODO: separate initialize==reset_state from __init__ 
 161  # TODO: introduce Ypos == diffC which makes the code more consistent and 
 162  #       the active update "exact"? 
 163  # TODO: dynamically read "signals" from a file, see myproperties.py 
 164  #               (to be called after tell()) 
 165  # 
 166  # typical parameters in scipy.optimize: disp, xtol, ftol, maxiter, maxfun, 
 167  #         callback=None 
 168  #         maxfev, diag (A sequency of N positive entries that serve as 
 169  #                 scale factors for the variables.) 
 170  #           full_output -- non-zero to return all optional outputs. 
 171  #   If xtol < 0.0, xtol is set to sqrt(machine_precision) 
 172  #    'infot -- a dictionary of optional outputs with the keys: 
 173  #                      'nfev': the number of function calls... 
 174  # 
 175  #    see eg fmin_powell 
 176  # typical returns 
 177  #        x, f, dictionary d 
 178  #        (xopt, {fopt, gopt, Hopt, func_calls, grad_calls, warnflag}, 
 179  #         <allvecs>) 
 180  # 
 181  # TODO: keep best ten solutions 
 182  # TODO: implement constraints handling 
 183  # TODO: extend function unitdoctest, or use unittest? 
 184  # TODO: apply style guide 
 185  # TODO: eigh(): thorough testing would not hurt 
 186   
 187  # changes: 
 188  # 14/05/07: added method result_pretty to pretty print optimization result 
 189  # 14/05/06: associated show() everywhere with ion() which should solve the 
 190  #           blocked terminal problem 
 191  # 14/05/05: all instances of "unicode" removed (was incompatible to 3.x) 
 192  # 14/05/05: replaced type(x) == y with isinstance(x, y), reorganized the 
 193  #           comments before the code starts 
 194  # 14/05/xx: change the order of kwargs of OOOptimizer.optimize, 
 195  #           remove prepare method in AdaptSigma classes, various changes/cleaning 
 196  # 14/03/01: bug fix BoundaryHandlerBase.has_bounds didn't check lower bounds correctly 
 197  #           bug fix in BoundPenalty.repair len(bounds[0]) was used instead of len(bounds[1]) 
 198  #           bug fix in GenoPheno.pheno, where x was not copied when only boundary-repair was applied 
 199  # 14/02/27: bug fixed when BoundPenalty was combined with fixed variables. 
 200  # 13/xx/xx: step-size adaptation becomes a class derived from CMAAdaptSigmaBase, 
 201  #           to make testing different adaptation rules (much) easier 
 202  # 12/12/14: separated CMAOptions and arguments to fmin 
 203  # 12/10/25: removed useless check_points from fmin interface 
 204  # 12/10/17: bug fix printing number of infeasible samples, moved not-in-use methods 
 205  #           timesCroot and divCroot to the right class 
 206  # 12/10/16 (0.92.00): various changes commit: bug bound[0] -> bounds[0], more_to_write fixed, 
 207  #   sigma_vec introduced, restart from elitist, trace normalization, max(mu,popsize/2) 
 208  #   is used for weight calculation. 
 209  # 12/07/23: (bug:) BoundPenalty.update respects now genotype-phenotype transformation 
 210  # 12/07/21: convert value True for noisehandling into 1 making the output compatible 
 211  # 12/01/30: class Solution and more old stuff removed r3101 
 212  # 12/01/29: class Solution is depreciated, GenoPheno and SolutionDict do the job (v0.91.00, r3100) 
 213  # 12/01/06: CMA_eigenmethod option now takes a function (integer still works) 
 214  # 11/09/30: flat fitness termination checks also history length 
 215  # 11/09/30: elitist option (using method clip_or_fit_solutions) 
 216  # 11/09/xx: method clip_or_fit_solutions for check_points option for all sorts of 
 217  #           injected or modified solutions and even reliable adaptive encoding 
 218  # 11/08/19: fixed: scaling and typical_x type clashes 1 vs array(1) vs ones(dim) vs dim * [1] 
 219  # 11/07/25: fixed: fmin wrote first and last line even with verb_log==0 
 220  #           fixed: method settableOptionsList, also renamed to versatileOptions 
 221  #           default seed depends on time now 
 222  # 11/07/xx (0.9.92): added: active CMA, selective mirrored sampling, noise/uncertainty handling 
 223  #           fixed: output argument ordering in fmin, print now only used as function 
 224  #           removed: parallel option in fmin 
 225  # 11/07/01: another try to get rid of the memory leak by replacing self.unrepaired = self[:] 
 226  # 11/07/01: major clean-up and reworking of abstract base classes and of the documentation, 
 227  #           also the return value of fmin changed and attribute stop is now a method. 
 228  # 11/04/22: bug-fix: option fixed_variables in combination with scaling 
 229  # 11/04/21: stopdict is not a copy anymore 
 230  # 11/04/15: option fixed_variables implemented 
 231  # 11/03/23: bug-fix boundary update was computed even without boundaries 
 232  # 11/03/12: bug-fix of variable annotation in plots 
 233  # 11/02/05: work around a memory leak in numpy 
 234  # 11/02/05: plotting routines improved 
 235  # 10/10/17: cleaning up, now version 0.9.30 
 236  # 10/10/17: bug-fix: return values of fmin now use phenotyp (relevant 
 237  #           if input scaling_of_variables is given) 
 238  # 08/10/01: option evalparallel introduced, 
 239  #           bug-fix for scaling being a vector 
 240  # 08/09/26: option CMAseparable becomes CMA_diagonal 
 241  # 08/10/18: some names change, test functions go into a class 
 242  # 08/10/24: more refactorizing 
 243  # 10/03/09: upper bound exp(min(1,...)) for step-size control 
 244   
 245  from __future__ import division 
 246  # future is >= 3.0, this code has mainly been used with 2.6 & 2.7 
 247  from __future__ import with_statement 
 248  # only necessary for python 2.5 and not in heavy use 
 249  from __future__ import print_function 
 250  # available from python 2.6, code should also work without 
 251  from __future__ import absolute_import 
 252  from __future__ import unicode_literals 
 253  # from __future__ import collections.MutableMapping 
 254  # does not exist in future, otherwise Python 2.5 would work, since 0.91.01 
 255   
 256  import sys 
 257  if not sys.version.startswith('2'):  # in python 3 
 258      xrange = range 
 259      raw_input = input 
 260      basestring = str 
 261  import time  # not really essential 
 262  import collections 
 263  import numpy as np 
 264  # arange, cos, size, eye, inf, dot, floor, outer, zeros, linalg.eigh, 
 265  # sort, argsort, random, ones,... 
 266  from numpy import inf, array, dot, exp, log, sqrt, sum 
 267  # to access the built-in sum fct:  ``__builtins__.sum`` or ``del sum`` 
 268  # removes the imported sum and recovers the shadowed build-in 
 269  try: 
 270      import matplotlib 
 271      import matplotlib.pyplot as pyplot  # also: use ipython -pyplot 
 272      savefig = pyplot.savefig  # now we can use cma.savefig() etc 
 273      closefig = pyplot.close 
 274 -    def show(): 
 275          # is_interactive = matplotlib.is_interactive() 
 276          pyplot.ion() 
 277          pyplot.show() 
 278          # if we call now matplotlib.interactive(True), the console is 
 279          # blocked 
 280      pyplot.ion()  # prevents that execution stops after plotting 
 281  except: 
 282      pyplot = None 
 283      savefig = None 
 284      closefig = None 
 285 -    def show(): 
 286          print('pyplot.show() is not available') 
 287      print('Could not import matplotlib.pyplot, therefore ``cma.plot()``" +' 
 288            ' etc. is not available') 
 289   
 290  __author__ = 'Nikolaus Hansen' 
 291  __version__ = "1.0.07  $Revision: 3809 $ $Date: 2014-05-08 02:37:39 +0200 (Thu, 08 May 2014) $" 
 292  # $Source$  # according to PEP 8 style guides, but what is it good for? 
 293  # $Id: cma.py 3809 2014-05-08 00:37:39Z hansen $ 
 294  # bash $: svn propset svn:keywords 'Date Revision Id' cma.py 
 295   
 296  __docformat__ = "reStructuredText"  # this hides some comments entirely? 
 297  __all__ = ( 
 298          'main', 
 299          'fmin', 
 300          'fcts', 
 301          'Fcts', 
 302          'felli', 
 303          'rotate', 
 304          'pprint', 
 305          'plot', 
 306          'disp', 
 307          'show', 
 308          'savefig', 
 309          'closefig', 
 310          'use_archives', 
 311          'is_feasible', 
 312          'unitdoctest', 
 313          'DerivedDictBase', 
 314          'SolutionDict', 
 315          'CMASolutionDict', 
 316          'BestSolution', 
 317          'BoundaryHandlerBase', 
 318          'BoundNone', 
 319          'BoundTransform', 
 320          'BoundPenalty', 
 321          'BoxConstraintsTransformationBase', 
 322          'BoxConstraintsLinQuadTransformation', 
 323          'GenoPheno', 
 324          'OOOptimizer', 
 325          'CMAEvolutionStrategy', 
 326          'CMAOptions', 
 327          'CMASolutionDict', 
 328          'CMAAdaptSigmaBase', 
 329          'CMAAdaptSigmaNone', 
 330          'CMAAdaptSigmaDistanceProportional', 
 331          'CMAAdaptSigmaCSA', 
 332          'CMAAdaptSigmaTPA', 
 333          'CMAAdaptSigmaMedianImprovement', 
 334          'BaseDataLogger', 
 335          'CMADataLogger', 
 336          'DEAPCMADataLogger', 
 337          'NoiseHandler', 
 338          'Sections', 
 339          'Misc', 
 340          'Mh', 
 341          'Rotation', 
 342          'FitnessFunctions' 
 343          ) 
 344  use_archives = True 
 345  # speed up for very large population size, prevents the need for an 
 346  # inverse gp-transformation, relies on collections module 
 347  # not sure what happens if set to False 
 348   
 349   
 350  # emptysets = ('', (), [], {}) 
 351  # array([]) does not work but np.size(.) == 0 
 352  # here is the problem: 
 353  # bool(array([0])) is False 
 354  # bool(list(array([0]))) is True 
 355  # bool(list(array([0, 1]))) is True 
 356  # bool(array([0, 1])) raises ValueError 
 357  # 
 358  # "x in emptysets" cannot be well replaced by "not x" 
 359  # which is also True for array([]) and None, but also for 0 and False, 
 360  # and False for NaN, and an exception for array([0,1]), see also 
 361  # http://google-styleguide.googlecode.com/svn/trunk/pyguide.html#True/False_evaluations 
 362   
 363  # ____________________________________________________________ 
 364  # ____________________________________________________________ 
 365  # 
 366 -def rglen(ar): 
 367      """shortcut for the iterator ``xrange(len(ar))``""" 
 368      return xrange(len(ar)) 
 369   
 370 -def is_feasible(x, f): 
 371      """default to check feasibility, see also ``cma_default_options``""" 
 372      return f is not None and f is not np.NaN 
 373   
 374 -def _print_warning(msg, method_name=None, class_name=None, iteration=None, 
 375                    verbose=1): 
 376      if verbose > 0: 
 377          print('WARNING (module=' + __name__ + 
 378                (', class=' + str(class_name) if class_name else '') + 
 379                (', method=' + str(method_name) if method_name else '') + 
 380                (', iteration=' + str(iteration) if iteration else '') + 
 381                '): ' + msg) 
 382   
 383  # ____________________________________________________________ 
 384  # ____________________________________________________________ 
 385  # 
 386 -def unitdoctest(): 
 387      """is used to describe test cases and might in future become helpful 
 388      as an experimental tutorial as well. The main testing feature at the 
 389      moment is by doctest with ``cma._test()`` or conveniently by 
 390      ``python cma.py --test``. With the ``--verbose`` option added, the 
 391      results will always slightly differ and many "failed" test cases 
 392      might be reported. 
 393   
 394      A simple first overall test: 
 395          >>> import cma 
 396          >>> res = cma.fmin(cma.fcts.elli, 3*[1], 1, 
 397          ...                {'CMA_diagonal':2, 'seed':1, 'verb_time':0}) 
 398          (3_w,7)-CMA-ES (mu_w=2.3,w_1=58%) in dimension 3 (seed=1) 
 399             Covariance matrix is diagonal for 2 iterations (1/ccov=7.0) 
 400          Iterat #Fevals   function value     axis ratio  sigma   minstd maxstd min:sec 
 401              1       7 1.453161670768570e+04 1.2e+00 1.08e+00  1e+00  1e+00 
 402              2      14 3.281197961927601e+04 1.3e+00 1.22e+00  1e+00  2e+00 
 403              3      21 1.082851071704020e+04 1.3e+00 1.24e+00  1e+00  2e+00 
 404            100     700 8.544042012075362e+00 1.4e+02 3.18e-01  1e-03  2e-01 
 405            200    1400 5.691152415221861e-12 1.0e+03 3.82e-05  1e-09  1e-06 
 406            220    1540 3.890107746209078e-15 9.5e+02 4.56e-06  8e-11  7e-08 
 407          termination on tolfun : 1e-11 
 408          final/bestever f-value = 3.89010774621e-15 2.52273602735e-15 
 409          mean solution:  [ -4.63614606e-08  -3.42761465e-10   1.59957987e-11] 
 410          std deviation: [  6.96066282e-08   2.28704425e-09   7.63875911e-11] 
 411   
 412      Test on the Rosenbrock function with 3 restarts. The first trial only 
 413      finds the local optimum, which happens in about 20% of the cases. 
 414   
 415          >>> import cma 
 416          >>> res = cma.fmin(cma.fcts.rosen, 4*[-1], 1, 
 417          ...                {'ftarget':1e-6, 'restarts':3, 
 418          ...                 'verb_time':0, 'verb_disp':500, 'seed':3}) 
 419          (4_w,8)-CMA-ES (mu_w=2.6,w_1=52%) in dimension 4 (seed=3) 
 420          Iterat #Fevals   function value     axis ratio  sigma   minstd maxstd min:sec 
 421              1       8 4.875315645656848e+01 1.0e+00 8.43e-01  8e-01  8e-01 
 422              2      16 1.662319948123120e+02 1.1e+00 7.67e-01  7e-01  8e-01 
 423              3      24 6.747063604799602e+01 1.2e+00 7.08e-01  6e-01  7e-01 
 424            184    1472 3.701428610430019e+00 4.3e+01 9.41e-07  3e-08  5e-08 
 425          termination on tolfun : 1e-11 
 426          final/bestever f-value = 3.70142861043 3.70142861043 
 427          mean solution:  [-0.77565922  0.61309336  0.38206284  0.14597202] 
 428          std deviation: [  2.54211502e-08   3.88803698e-08   4.74481641e-08   3.64398108e-08] 
 429          (8_w,16)-CMA-ES (mu_w=4.8,w_1=32%) in dimension 4 (seed=4) 
 430          Iterat #Fevals   function value     axis ratio  sigma   minstd maxstd min:sec 
 431              1    1489 2.011376859371495e+02 1.0e+00 8.90e-01  8e-01  9e-01 
 432              2    1505 4.157106647905128e+01 1.1e+00 8.02e-01  7e-01  7e-01 
 433              3    1521 3.548184889359060e+01 1.1e+00 1.02e+00  8e-01  1e+00 
 434            111    3249 6.831867555502181e-07 5.1e+01 2.62e-02  2e-04  2e-03 
 435          termination on ftarget : 1e-06 
 436          final/bestever f-value = 6.8318675555e-07 1.18576673231e-07 
 437          mean solution:  [ 0.99997004  0.99993938  0.99984868  0.99969505] 
 438          std deviation: [ 0.00018973  0.00038006  0.00076479  0.00151402] 
 439          >>> assert res[1] <= 1e-6 
 440   
 441      Notice the different termination conditions. Termination on the target 
 442      function value ftarget prevents further restarts. 
 443   
 444      Test of scaling_of_variables option 
 445   
 446          >>> import cma 
 447          >>> opts = cma.CMAOptions() 
 448          >>> opts['seed'] = 456 
 449          >>> opts['verb_disp'] = 0 
 450          >>> opts['CMA_active'] = 1 
 451          >>> # rescaling of third variable: for searching in  roughly 
 452          >>> #   x0 plus/minus 1e3*sigma0 (instead of plus/minus sigma0) 
 453          >>> opts.scaling_of_variables = [1, 1, 1e3, 1] 
 454          >>> res = cma.fmin(cma.fcts.rosen, 4 * [0.1], 0.1, opts) 
 455          termination on tolfun : 1e-11 
 456          final/bestever f-value = 2.68096173031e-14 1.09714829146e-14 
 457          mean solution:  [ 1.00000001  1.00000002  1.00000004  1.00000007] 
 458          std deviation: [  3.00466854e-08   5.88400826e-08   1.18482371e-07   2.34837383e-07] 
 459   
 460      The printed std deviations reflect the actual true value (not the one 
 461      in the internal representation which would be different). 
 462   
 463          >>> import cma 
 464          >>> r = cma.fmin(cma.fcts.diffpow, 15 * [1], 1, 
 465          ...              {'CMA_dampsvec_fac':0.5, 'ftarget':1e-9}) 
 466          >>> assert(r[1] < 1e-9) 
 467          >>> assert(r[2] < 13000)  # only passed with CMA_dampsvec_fac 
 468   
 469   
 470      :See: cma.main(), cma._test() 
 471   
 472      """ 
 473   
 474      pass 
 475   
 476  # ____________________________________________________________ 
 477  # ____________________________________________________________ 
 478  # 
 479 -class _BlancClass(object): 
 480      """blanc container class for having a collection of attributes""" 
 481   
 482  # _____________________________________________________________________ 
 483  # _____________________________________________________________________ 
 484  # 
 485 -class DerivedDictBase(collections.MutableMapping): 
 486      """for conveniently adding "features" to a dictionary. The actual 
 487      dictionary is in ``self.data``. Copy-paste 
 488      and modify setitem, getitem, and delitem, if necessary""" 
 489 -    def __init__(self, *args, **kwargs): 
 490          # collections.MutableMapping.__init__(self) 
 491          super(DerivedDictBase, self).__init__() 
 492          # super(SolutionDict, self).__init__()  # the same 
 493          self.data = dict(*args, **kwargs) 
 494 -    def __len__(self): 
 495          return len(self.data) 
 496 -    def __contains__(self, value): 
 497          return value in self.data 
 498 -    def __iter__(self): 
 499          return iter(self.data) 
 500 -    def __setitem__(self, key, value): 
 501          """defines self[key] = value""" 
 502          self.data[key] = value 
 503 -    def __getitem__(self, key): 
 504          """defines self[key]""" 
 505          return self.data[key] 
 506 -    def __delitem__(self, key): 
 507          del self.data[key] 
 508   
 509 -class SolutionDict(DerivedDictBase): 
 510      """dictionary with computation of an hash key. 
 511   
 512      The hash key is generated from the inserted solution and a stack of 
 513      previously inserted same solutions is provided. Each entry is meant 
 514      to store additional information related to the solution. 
 515   
 516          >>> import cma, numpy as np 
 517          >>> d = cma.SolutionDict() 
 518          >>> x = np.array([1,2,4]) 
 519          >>> d[x] = {'f': sum(x**2), 'iteration': 1} 
 520          >>> assert d[x]['iteration'] == 1 
 521          >>> assert d.get(x) == (d[x] if d.key(x) in d.keys() else None) 
 522   
 523      TODO: data_with_same_key behaves like a stack (see setitem and 
 524      delitem), but rather should behave like a queue?! A queue is less 
 525      consistent with the operation self[key] = ..., if 
 526      self.data_with_same_key[key] is not empty. 
 527   
 528      TODO: iteration key is used to clean up without error management 
 529   
 530      """ 
 531 -    def __init__(self, *args, **kwargs): 
 532          DerivedDictBase.__init__(self, *args, **kwargs) 
 533          self.data_with_same_key = {} 
 534          self.last_iteration = 0 
 535 -    def key(self, x): 
 536          try: 
 537              return tuple(x) 
 538              # using sum(x) is slower, using x[0] is slightly faster 
 539          except TypeError: 
 540              return x 
 541 -    def __setitem__(self, key, value): 
 542          """defines self[key] = value""" 
 543          key = self.key(key) 
 544          if key in self.data_with_same_key: 
 545              self.data_with_same_key[key] += [self.data[key]] 
 546          elif key in self.data: 
 547              self.data_with_same_key[key] = [self.data[key]] 
 548          self.data[key] = value 
 549 -    def __getitem__(self, key):  # 50% of time of 
 550          """defines self[key]""" 
 551          return self.data[self.key(key)] 
 552 -    def __delitem__(self, key): 
 553          """remove only most current key-entry""" 
 554          key = self.key(key) 
 555          if key in self.data_with_same_key: 
 556              if len(self.data_with_same_key[key]) == 1: 
 557                  self.data[key] = self.data_with_same_key.pop(key)[0] 
 558              else: 
 559                  self.data[key] = self.data_with_same_key[key].pop(-1) 
 560          else: 
 561              del self.data[key] 
 562 -    def truncate(self, max_len, min_iter): 
 563          if len(self) > max_len: 
 564              for k in list(self.keys()): 
 565                  if self[k]['iteration'] < min_iter: 
 566                      del self[k] 
 567                      # deletes one item with k as key, better delete all? 
 568   
 569 -class CMASolutionDict(SolutionDict): 
 570 -    def __init__(self, *args, **kwargs): 
 571          SolutionDict.__init__(self, *args, **kwargs) 
 572          self.last_solution_index = 0 
 573   
 574      # TODO: insert takes 30% of the overall CPU time, mostly in def key() 
 575      #       with about 15% of the overall CPU time 
 576 -    def insert(self, key, geno=None, iteration=None, fitness=None, value=None): 
 577          """insert an entry with key ``key`` and value 
 578          ``value if value is not None else {'geno':key}`` and 
 579          ``self[key]['kwarg'] = kwarg if kwarg is not None`` for the further kwargs. 
 580   
 581          """ 
 582          # archive returned solutions, first clean up archive 
 583          if iteration is not None and iteration > self.last_iteration and (iteration % 10) < 1: 
 584              self.truncate(300, iteration - 3) 
 585          elif value is not None and value.get('iteration'): 
 586              iteration = value['iteration'] 
 587              if (iteration % 10) < 1: 
 588                  self.truncate(300, iteration - 3) 
 589   
 590          self.last_solution_index += 1 
 591          if value is not None: 
 592              try: 
 593                  iteration = value['iteration'] 
 594              except: 
 595                  pass 
 596          if iteration is not None: 
 597              if iteration > self.last_iteration: 
 598                  self.last_solution_index = 0 
 599              self.last_iteration = iteration 
 600          else: 
 601              iteration = self.last_iteration + 0.5  # a hack to get a somewhat reasonable value 
 602          if value is not None: 
 603              self[key] = value 
 604          else: 
 605              self[key] = {'pheno': key} 
 606          if geno is not None: 
 607              self[key]['geno'] = geno 
 608          if iteration is not None: 
 609              self[key]['iteration'] = iteration 
 610          if fitness is not None: 
 611              self[key]['fitness'] = fitness 
 612          return self[key] 
 613   
 614  if not use_archives: 
 615 -    class CMASolutionDict(SolutionDict): 
 616 -        def insert(self, *args, **kwargs): 
 617              pass 
 618   
 619 -class BestSolution(object): 
 620      """container to keep track of the best solution seen""" 
 621 -    def __init__(self, x=None, f=np.inf, evals=None): 
 622          """initialize the best solution with `x`, `f`, and `evals`. 
 623          Better solutions have smaller `f`-values. 
 624   
 625          """ 
 626          self.x = x 
 627          self.x_geno = None 
 628          self.f = f if f is not None and f is not np.nan else np.inf 
 629          self.evals = evals 
 630          self.evalsall = evals 
 631          self.last = _BlancClass() 
 632          self.last.x = x 
 633          self.last.f = f 
 634 -    def update(self, arx, xarchive=None, arf=None, evals=None): 
 635          """checks for better solutions in list `arx`. 
 636   
 637          Based on the smallest corresponding value in `arf`, 
 638          alternatively, `update` may be called with a `BestSolution` 
 639          instance like ``update(another_best_solution)`` in which case 
 640          the better solution becomes the current best. 
 641   
 642          `xarchive` is used to retrieve the genotype of a solution. 
 643   
 644          """ 
 645          if isinstance(arx, BestSolution): 
 646              if self.evalsall is None: 
 647                  self.evalsall = arx.evalsall 
 648              elif arx.evalsall is not None: 
 649                  self.evalsall = max((self.evalsall, arx.evalsall)) 
 650              if arx.f is not None and arx.f < np.inf: 
 651                  self.update([arx.x], xarchive, [arx.f], arx.evals) 
 652              return self 
 653          assert arf is not None 
 654          # find failsave minimum 
 655          minidx = np.nanargmin(arf) 
 656          if minidx is np.nan: 
 657              return 
 658          minarf = arf[minidx] 
 659          # minarf = reduce(lambda x, y: y if y and y is not np.nan 
 660          #                   and y < x else x, arf, np.inf) 
 661          if minarf < np.inf and (minarf < self.f or self.f is None): 
 662              self.x, self.f = arx[minidx], arf[minidx] 
 663              if xarchive is not None and xarchive.get(self.x) is not None: 
 664                  self.x_geno = xarchive[self.x].get('geno') 
 665              else: 
 666                  self.x_geno = None 
 667              self.evals = None if not evals else evals - len(arf) + minidx + 1 
 668              self.evalsall = evals 
 669          elif evals: 
 670              self.evalsall = evals 
 671          self.last.x = arx[minidx] 
 672          self.last.f = minarf 
 673 -    def get(self): 
 674          """return ``(x, f, evals)`` """ 
 675          return self.x, self.f, self.evals  # , self.x_geno 
 676   
 677   
 678  # ____________________________________________________________ 
 679  # ____________________________________________________________ 
 680  # 
 681 -class BoundaryHandlerBase(object): 
 682      """hacked base class """ 
 683 -    def __init__(self, bounds): 
 684          """bounds are not copied, but possibly modified and 
 685          put into a normalized form: ``bounds`` can be ``None`` 
 686          or ``[lb, ub]`` where ``lb`` and ``ub`` are 
 687          either None or a vector (which can have ``None`` entries). 
 688   
 689          Generally, the last entry is recycled to compute bounds 
 690          for any dimension. 
 691   
 692          """ 
 693          if not bounds: 
 694              self.bounds = None 
 695          else: 
 696              l = [None, None]  # figure out lenths 
 697              for i in [0, 1]: 
 698                  try: 
 699                      l[i] = len(bounds[i]) 
 700                  except TypeError: 
 701                      bounds[i] = [bounds[i]] 
 702                      l[i] = 1 
 703                  if all([bounds[i][j] is None or not np.isfinite(bounds[i][j]) 
 704                          for j in rglen(bounds[i])]): 
 705                      bounds[i] = None 
 706                  if bounds[i] is not None and any([bounds[i][j] == (-1)**i * np.inf 
 707                                                    for j in rglen(bounds[i])]): 
 708                      raise ValueError('lower/upper is +inf/-inf and ' + 
 709                                       'therefore no finite feasible solution is available') 
 710              self.bounds = bounds 
 711   
 712 -    def __call__(self, solutions, *args, **kwargs): 
 713          """return penalty or list of penalties, by default zero(s). 
 714   
 715          This interface seems too specifically tailored to the derived 
 716          BoundPenalty class, it should maybe change. 
 717   
 718          """ 
 719          if np.isscalar(solutions[0]): 
 720              return 0.0 
 721          else: 
 722              return len(solutions) * [0.0] 
 723   
 724 -    def update(self, *args, **kwargs): 
 725          return self 
 726   
 727 -    def repair(self, x, copy_if_changed=True, copy_always=False): 
 728          """projects infeasible values on the domain bound, might be 
 729          overwritten by derived class """ 
 730          if copy_always: 
 731              x = array(x, copy=True) 
 732              copy = False 
 733          else: 
 734              copy = copy_if_changed 
 735          if self.bounds is None: 
 736              return x 
 737          for ib in [0, 1]: 
 738              if self.bounds[ib] is None: 
 739                  continue 
 740              for i in rglen(x): 
 741                  idx = min([i, len(self.bounds[ib]) - 1]) 
 742                  if self.bounds[ib][idx] is not None and \ 
 743                          (-1)**ib * x[i] < (-1)**ib * self.bounds[ib][idx]: 
 744                      if copy: 
 745                          x = array(x, copy=True) 
 746                          copy = False 
 747                      x[i] = self.bounds[ib][idx] 
 748   
 749 -    def inverse(self, y, copy_if_changed=True, copy_always=False): 
 750          return y if not copy_always else array(y, copy=True) 
 751   
 752 -    def get_bounds(self, which, dimension): 
 753          """``get_bounds('lower', 8)`` returns the lower bounds in 8-D""" 
 754          if which == 'lower' or which == 0: 
 755              return self._get_bounds(0, dimension) 
 756          elif which == 'upper' or which == 1: 
 757              return self._get_bounds(1, dimension) 
 758          else: 
 759              raise ValueError("argument which must be 'lower' or 'upper'") 
 760   
 761 -    def _get_bounds(self, ib, dimension): 
 762          """ib == 0/1 means lower/upper bound, return a vector of length 
 763          `dimension` """ 
 764          sign_ = 2 * ib - 1 
 765          assert sign_**2 == 1 
 766          if self.bounds is None or self.bounds[ib] is None: 
 767              return array(dimension * [sign_ * np.Inf]) 
 768          res = [] 
 769          for i in xrange(dimension): 
 770              res.append(self.bounds[ib][min([i, len(self.bounds[ib]) - 1])]) 
 771              if res[-1] is None: 
 772                  res[-1] = sign_ * np.Inf 
 773          return array(res) 
 774   
 775 -    def has_bounds(self): 
 776          """return True, if any variable is bounded""" 
 777          bounds = self.bounds 
 778          if bounds in (None, [None, None]): 
 779              return False 
 780          for ib, bound in enumerate(bounds): 
 781              if bound is not None: 
 782                  sign_ = 2 * ib - 1 
 783                  for bound_i in bound: 
 784                      if bound_i is not None and sign_ * bound_i < np.inf: 
 785                          return True 
 786          return False 
 787   
 788 -    def is_in_bounds(self, x): 
 789          """not yet tested""" 
 790          if self.bounds is None: 
 791              return True 
 792          for ib in [0, 1]: 
 793              if self.bounds[ib] is None: 
 794                  continue 
 795              for i in rglen(x): 
 796                  idx = min([i, len(self.bounds[ib]) - 1]) 
 797                  if self.bounds[ib][idx] is not None and \ 
 798                          (-1)**ib * x[i] < (-1)**ib * self.bounds[ib][idx]: 
 799                      return False 
 800          return True 
 801   
 802 -    def to_dim_times_two(self, bounds): 
 803          """return boundaries in format ``[[lb0, ub0], [lb1, ub1], ...]``, 
 804          as used by ``BoxConstraints...`` class. 
 805   
 806          """ 
 807          if not bounds: 
 808              b = [[None, None]] 
 809          else: 
 810              l = [None, None]  # figure out lenths 
 811              for i in [0, 1]: 
 812                  try: 
 813                      l[i] = len(bounds[i]) 
 814                  except TypeError: 
 815                      bounds[i] = [bounds[i]] 
 816                      l[i] = 1 
 817              b = []  # bounds in different format 
 818              try: 
 819                  for i in xrange(max(l)): 
 820                      b.append([bounds[0][i] if i < l[0] else None, 
 821                                bounds[1][i] if i < l[1] else None]) 
 822              except (TypeError, IndexError): 
 823                  print("boundaries must be provided in the form " + 
 824                        "[scalar_of_vector, scalar_or_vector]") 
 825                  raise 
 826          return b 
 827   
 828  # ____________________________________________________________ 
 829  # ____________________________________________________________ 
 830  # 
 831 -class BoundNone(BoundaryHandlerBase): 
 832 -    def __init__(self, bounds=None): 
 833          if bounds is not None: 
 834              raise ValueError() 
 835          BoundaryHandlerBase.__init__(self, None) 
 836 -    def is_in_bounds(self, x): 
 837          return True 
 838   
 839  # ____________________________________________________________ 
 840  # ____________________________________________________________ 
 841  # 
 842 -class BoundTransform(BoundaryHandlerBase): 
 843      """Handles boundary by a smooth, piecewise linear and quadratic 
 844      transformation into the feasible domain. 
 845   
 846      >>> import cma 
 847      >>> b = cma.BoundTransform([None, 1]) 
 848      >>> assert b.bounds == [[None], [1]] 
 849      >>> assert cma.Mh.vequals_approximately(b.repair([0, 1, 1.2]), 
 850      ...                                     array([ 0., 0.975, 0.975])) 
 851      >>> assert b.is_in_bounds([0, 0.5, 1]) 
 852      >>> assert cma.Mh.vequals_approximately(b.transform([0, 1, 2]), 
 853      ...                                     [ 0.   ,  0.975,  0.2  ]) 
 854      >>> o=cma.fmin(cma.fcts.sphere, 6 * [-2], 0.5, options={ 
 855      ...    'boundary_handling': 'BoundTransform ', 
 856      ...    'bounds': [[], 5 * [-1] + [inf]] }) 
 857      >>> assert o[1] < 5 + 1e-8 
 858   
 859      Details: this class uses ``class BoxConstraintsLinQuadTransformation`` 
 860   
 861      """ 
 862 -    def __init__(self, bounds=None): 
 863          """Argument bounds can be `None` or ``bounds[0]`` and ``bounds[1]`` 
 864          are lower and upper domain boundaries, each is either `None` or 
 865          a scalar or a list or array of appropriate size. 
 866   
 867          """ 
 868          BoundaryHandlerBase.__init__(self, bounds) 
 869          self.bounds_tf = BoxConstraintsLinQuadTransformation(self.to_dim_times_two(bounds)) 
 870   
 871 -    def repair(self, x, copy_if_changed=True, copy_always=False): 
 872          """transforms ``x`` into the bounded domain. 
 873   
 874          ``copy_always`` option might disappear. 
 875   
 876          """ 
 877          copy = copy_if_changed 
 878          if copy_always: 
 879              x = array(x, copy=True) 
 880              copy = False 
 881          if self.bounds is None or (self.bounds[0] is None and 
 882                                     self.bounds[1] is None): 
 883              return x 
 884          return self.bounds_tf(x, copy) 
 885   
 886 -    def transform(self, x): 
 887          return self.repair(x) 
 888   
 889 -    def inverse(self, x, copy_if_changed=True, copy_always=False): 
 890          """inverse transform of ``x`` from the bounded domain. 
 891   
 892          """ 
 893          copy = copy_if_changed 
 894          if copy_always: 
 895              x = array(x, copy=True) 
 896              copy = False 
 897          if self.bounds is None or (self.bounds[0] is None and 
 898                                     self.bounds[1] is None): 
 899              return x 
 900          return self.bounds_tf.inverse(x, copy)  # this doesn't exist 
 901   
 902  # ____________________________________________________________ 
 903  # ____________________________________________________________ 
 904  # 
 905 -class BoundPenalty(BoundaryHandlerBase): 
 906      """Computes the boundary penalty. Must be updated each iteration, 
 907      using the `update` method. 
 908   
 909      Details 
 910      ------- 
 911      The penalty computes like ``sum(w[i] * (x[i]-xfeas[i])**2)``, 
 912      where `xfeas` is the closest feasible (in-bounds) solution from `x`. 
 913      The weight `w[i]` should be updated during each iteration using 
 914      the update method. 
 915   
 916      Example: 
 917   
 918      >>> import cma 
 919      >>> cma.fmin(cma.felli, 6 * [1], 1, 
 920      ...          { 
 921      ...              'boundary_handling': 'BoundPenalty', 
 922      ...              'bounds': [-1, 1], 
 923      ...              'fixed_variables': {0: 0.012, 2:0.234} 
 924      ...          }) 
 925   
 926      Reference: Hansen et al 2009, A Method for Handling Uncertainty... 
 927      IEEE TEC, with addendum, see 
 928      http://www.lri.fr/~hansen/TEC2009online.pdf 
 929   
 930      """ 
 931 -    def __init__(self, bounds=None): 
 932          """Argument bounds can be `None` or ``bounds[0]`` and ``bounds[1]`` 
 933          are lower  and upper domain boundaries, each is either `None` or 
 934          a scalar or a list or array of appropriate size. 
 935          """ 
 936          # # 
 937          # bounds attribute reminds the domain boundary values 
 938          BoundaryHandlerBase.__init__(self, bounds) 
 939   
 940          self.gamma = 1  # a very crude assumption 
 941          self.weights_initialized = False  # gamma becomes a vector after initialization 
 942          self.hist = []  # delta-f history 
 943   
 944 -    def repair(self, x, copy_if_changed=True, copy_always=False): 
 945          """sets out-of-bounds components of ``x`` on the bounds. 
 946   
 947          """ 
 948          # TODO (old data): CPU(N,lam,iter=20,200,100): 3.3s of 8s for two bounds, 1.8s of 6.5s for one bound 
 949          # remark: np.max([bounds[0], x]) is about 40 times slower than max((bounds[0], x)) 
 950          copy = copy_if_changed 
 951          if copy_always: 
 952              x = array(x, copy=True) 
 953          bounds = self.bounds 
 954          if bounds not in (None, [None, None], (None, None)):  # solely for effiency 
 955              x = array(x, copy=True) if copy and not copy_always else x 
 956              if bounds[0] is not None: 
 957                  if np.isscalar(bounds[0]): 
 958                      for i in rglen(x): 
 959                          x[i] = max((bounds[0], x[i])) 
 960                  else: 
 961                      for i in rglen(x): 
 962                          j = min([i, len(bounds[0]) - 1]) 
 963                          if bounds[0][j] is not None: 
 964                              x[i] = max((bounds[0][j], x[i])) 
 965              if bounds[1] is not None: 
 966                  if np.isscalar(bounds[1]): 
 967                      for i in rglen(x): 
 968                          x[i] = min((bounds[1], x[i])) 
 969                  else: 
 970                      for i in rglen(x): 
 971                          j = min((i, len(bounds[1]) - 1)) 
 972                          if bounds[1][j] is not None: 
 973                              x[i] = min((bounds[1][j], x[i])) 
 974          return x 
 975   
 976      # ____________________________________________________________ 
 977      # 
 978 -    def __call__(self, x, archive, gp): 
 979          """returns the boundary violation penalty for `x` ,where `x` is a 
 980          single solution or a list or array of solutions. 
 981   
 982          """ 
 983          if x in (None, (), []): 
 984              return x 
 985          if self.bounds in (None, [None, None], (None, None)): 
 986              return 0.0 if np.isscalar(x[0]) else [0.0] * len(x)  # no penalty 
 987   
 988          x_is_single_vector = np.isscalar(x[0]) 
 989          x = [x] if x_is_single_vector else x 
 990   
 991          # add fixed variables to self.gamma 
 992          try: 
 993              gamma = list(self.gamma)  # fails if self.gamma is a scalar 
 994              for i in sorted(gp.fixed_values):  # fails if fixed_values is None 
 995                  gamma.insert(i, 0.0) 
 996              gamma = array(gamma, copy=False) 
 997          except TypeError: 
 998              gamma = self.gamma 
 999          pen = [] 
1000          for xi in x: 
1001              # CAVE: this does not work with already repaired values!! 
1002              # CPU(N,lam,iter=20,200,100)?: 3s of 10s, array(xi): 1s 
1003              # remark: one deep copy can be prevented by xold = xi first 
1004              xpheno = gp.pheno(archive[xi]['geno']) 
1005              # necessary, because xi was repaired to be in bounds 
1006              xinbounds = self.repair(xpheno) 
1007              # could be omitted (with unpredictable effect in case of external repair) 
1008              fac = 1  # exp(0.1 * (log(self.scal) - np.mean(self.scal))) 
1009              pen.append(sum(gamma * ((xinbounds - xpheno) / fac)**2) / len(xi)) 
1010          return pen[0] if x_is_single_vector else pen 
1011   
1012      # ____________________________________________________________ 
1013      # 
1014 -    def feasible_ratio(self, solutions): 
1015          """counts for each coordinate the number of feasible values in 
1016          ``solutions`` and returns an array of length ``len(solutions[0])`` 
1017          with the ratios. 
1018   
1019          `solutions` is a list or array of repaired `Solution` instances 
1020   
1021          """ 
1022          count = np.zeros(len(solutions[0])) 
1023          for x in solutions: 
1024              count += x.unrepaired == x 
1025          return count / float(len(solutions)) 
1026   
1027      # ____________________________________________________________ 
1028      # 
1029 -    def update(self, function_values, es): 
1030          """updates the weights for computing a boundary penalty. 
1031   
1032          Arguments 
1033          --------- 
1034          `function_values` 
1035              all function values of recent population of solutions 
1036          `es` 
1037              `CMAEvolutionStrategy` object instance, in particular 
1038              mean and variances and the methods from the attribute 
1039              `gp` of type `GenoPheno` are used. 
1040   
1041          """ 
1042          if self.bounds is None or (self.bounds[0] is None and 
1043                                     self.bounds[1] is None): 
1044              return self 
1045   
1046          N = es.N 
1047          # ## prepare 
1048          # compute varis = sigma**2 * C_ii 
1049          varis = es.sigma**2 * array(N * [es.C] if np.isscalar(es.C) else (# scalar case 
1050                                  es.C if np.isscalar(es.C[0]) else  # diagonal matrix case 
1051                                  [es.C[i][i] for i in xrange(N)]))  # full matrix case 
1052   
1053          # relative violation in geno-space 
1054          dmean = (es.mean - es.gp.geno(self.repair(es.gp.pheno(es.mean)))) / varis**0.5 
1055   
1056          # ## Store/update a history of delta fitness value 
1057          fvals = sorted(function_values) 
1058          l = 1 + len(fvals) 
1059          val = fvals[3 * l // 4] - fvals[l // 4]  # exact interquartile range apart interpolation 
1060          val = val / np.mean(varis)  # new: val is normalized with sigma of the same iteration 
1061          # insert val in history 
1062          if np.isfinite(val) and val > 0: 
1063              self.hist.insert(0, val) 
1064          elif val == inf and len(self.hist) > 1: 
1065              self.hist.insert(0, max(self.hist)) 
1066          else: 
1067              pass  # ignore 0 or nan values 
1068          if len(self.hist) > 20 + (3 * N) / es.popsize: 
1069              self.hist.pop() 
1070   
1071          # ## prepare 
1072          dfit = np.median(self.hist)  # median interquartile range 
1073          damp = min(1, es.sp.mueff / 10. / N) 
1074   
1075          # ## set/update weights 
1076          # Throw initialization error 
1077          if len(self.hist) == 0: 
1078              raise _Error('wrongful initialization, no feasible solution sampled. ' + 
1079                  'Reasons can be mistakenly set bounds (lower bound not smaller than upper bound) or a too large initial sigma0 or... ' + 
1080                  'See description of argument func in help(cma.fmin) or an example handling infeasible solutions in help(cma.CMAEvolutionStrategy). ') 
1081          # initialize weights 
1082          if (dmean.any() and (not self.weights_initialized or es.countiter == 2)):  # TODO 
1083              self.gamma = array(N * [2 * dfit])  ## BUGBUGzzzz: N should be phenotypic (bounds are in phenotype), but is genotypic 
1084              self.weights_initialized = True 
1085          # update weights gamma 
1086          if self.weights_initialized: 
1087              edist = array(abs(dmean) - 3 * max(1, N**0.5 / es.sp.mueff)) 
1088              if 1 < 3:  # this is better, around a factor of two 
1089                  # increase single weights possibly with a faster rate than they can decrease 
1090                  #     value unit of edst is std dev, 3==random walk of 9 steps 
1091                  self.gamma *= exp((edist > 0) * np.tanh(edist / 3) / 2.)**damp 
1092                  # decrease all weights up to the same level to avoid single extremely small weights 
1093                  #    use a constant factor for pseudo-keeping invariance 
1094                  self.gamma[self.gamma > 5 * dfit] *= exp(-1. / 3)**damp 
1095                  #     self.gamma[idx] *= exp(5*dfit/self.gamma[idx] - 1)**(damp/3) 
1096          es.more_to_write += list(self.gamma) if self.weights_initialized else N * [1.0] 
1097          # ## return penalty 
1098          # es.more_to_write = self.gamma if not np.isscalar(self.gamma) else N*[1] 
1099          return self  # bound penalty values 
1100   
1101  # ____________________________________________________________ 
1102  # ____________________________________________________________ 
1103  # 
1104 -class BoxConstraintsTransformationBase(object): 
1105      """Implements a transformation into boundaries and is used for 
1106      boundary handling:: 
1107   
1108          tf = BoxConstraintsTransformationAnyDerivedClass([[1, 4]]) 
1109          x = [3, 2, 4.4] 
1110          y = tf(x)  # "repaired" solution 
1111          print(tf([2.5]))  # middle value is never changed 
1112          [2.5] 
1113   
1114      :See: ``BoundaryHandler`` 
1115   
1116      """ 
1117 -    def __init__(self, bounds): 
1118          try: 
1119              if len(bounds[0]) != 2: 
1120                  raise ValueError 
1121          except: 
1122              raise ValueError(' bounds must be either [[lb0, ub0]] or [[lb0, ub0], [lb1, ub1],...], \n where in both cases the last entry is reused for all remaining dimensions') 
1123          self.bounds = bounds 
1124          self.initialize() 
1125   
1126 -    def initialize(self): 
1127          """initialize in base class""" 
1128          self._lb = [b[0] for b in self.bounds]  # can be done more efficiently? 
1129          self._ub = [b[1] for b in self.bounds] 
1130   
1131 -    def _lowerupperval(self, a, b, c): 
1132          return np.max([np.max(a), np.min([np.min(b), c])]) 
1133 -    def bounds_i(self, i): 
1134          """return ``[ith_lower_bound, ith_upper_bound]``""" 
1135          return self.bounds[self._index(i)] 
1136 -    def __call__(self, solution_in_genotype): 
1137          res = [self._transform_i(x, i) for i, x in enumerate(solution_in_genotype)] 
1138          return res 
1139      transform = __call__ 
1140 -    def inverse(self, solution_in_phenotype, copy_if_changed=True, copy_always=True): 
1141          return [self._inverse_i(y, i) for i, y in enumerate(solution_in_phenotype)] 
1142 -    def _index(self, i): 
1143          return min((i, len(self.bounds) - 1)) 
1144 -    def _transform_i(self, x, i): 
1145          raise NotImplementedError('this is an abstract method that should be implemented in the derived class') 
1146 -    def _inverse_i(self, y, i): 
1147          raise NotImplementedError('this is an abstract method that should be implemented in the derived class') 
1148 -    def shift_or_mirror_into_invertible_domain(self, solution_genotype): 
1149          """return the reference solution that has the same ``box_constraints_transformation(solution)`` 
1150          value, i.e. ``tf.shift_or_mirror_into_invertible_domain(x) = tf.inverse(tf.transform(x))``. 
1151          This is an idempotent mapping (leading to the same result independent how often it is 
1152          repeatedly applied). 
1153   
1154          """ 
1155          return self.inverse(self(solution_genotype)) 
1156          raise NotImplementedError('this is an abstract method that should be implemented in the derived class') 
1157   
1158 -class _BoxConstraintsTransformationTemplate(BoxConstraintsTransformationBase): 
1159      """copy/paste this template to implement a new boundary handling transformation""" 
1160 -    def __init__(self, bounds): 
1161          BoxConstraintsTransformationBase.__init__(self, bounds) 
1162 -    def initialize(self): 
1163          BoxConstraintsTransformationBase.initialize(self)  # likely to be removed 
1164 -    def _transform_i(self, x, i): 
1165          raise NotImplementedError('this is an abstract method that should be implemented in the derived class') 
1166 -    def _inverse_i(self, y, i): 
1167          raise NotImplementedError('this is an abstract method that should be implemented in the derived class') 
1168      __doc__ = BoxConstraintsTransformationBase.__doc__ + __doc__ 
1169   
1170 -class BoxConstraintsLinQuadTransformation(BoxConstraintsTransformationBase): 
1171      """implements a bijective, monotonous transformation between [lb - al, ub + au] 
1172      and [lb, ub] which is the identity (and therefore linear) in [lb + al, ub - au] 
1173      (typically about 90% of the interval) and quadratic in [lb - 3*al, lb + al] 
1174      and in [ub - au, ub + 3*au]. The transformation is periodically 
1175      expanded beyond the limits (somewhat resembling the shape sin(x-pi/2)) 
1176      with a period of ``2 * (ub - lb + al + au)``. 
1177   
1178      Details 
1179      ======= 
1180      Partly due to numerical considerations depend the values ``al`` and ``au`` 
1181      on ``abs(lb)`` and ``abs(ub)`` which makes the transformation non-translation 
1182      invariant. In contrast to sin(.), the transformation is robust to "arbitrary" 
1183      values for boundaries, e.g. a lower bound of ``-1e99`` or ``np.Inf`` or 
1184      ``None``. 
1185   
1186      Examples 
1187      ======== 
1188      Example to use with cma: 
1189   
1190      >>> import cma 
1191      >>> # only the first variable has an upper bound 
1192      >>> tf = cma.BoxConstraintsLinQuadTransformation([[1,2], [1,None]]) # second==last pair is re-cycled 
1193      >>> cma.fmin(cma.felli, 9 * [2], 1, {'transformation': [tf.transform, tf.inverse], 'verb_disp': 0}) 
1194      >>> # ...or... 
1195      >>> es = cma.CMAEvolutionStrategy(9 * [2], 1) 
1196      >>> while not es.stop(): 
1197      ...     X = es.ask() 
1198      ...     f = [cma.felli(tf(x)) for x in X]  # tf(x) == tf.transform(x) 
1199      ...     es.tell(X, f) 
1200   
1201      Example of the internal workings: 
1202   
1203      >>> import cma 
1204      >>> tf = cma.BoxConstraintsLinQuadTransformation([[1,2], [1,11], [1,11]]) 
1205      >>> tf.bounds 
1206      [[1, 2], [1, 11], [1, 11]] 
1207      >>> tf([1.5, 1.5, 1.5]) 
1208      [1.5, 1.5, 1.5] 
1209      >>> tf([1.52, -2.2, -0.2, 2, 4, 10.4]) 
1210      [1.52, 4.0, 2.0, 2.0, 4.0, 10.4] 
1211      >>> res = np.round(tf._au, 2) 
1212      >>> assert list(res[:4]) == [ 0.15, 0.6, 0.6, 0.6] 
1213      >>> res = [round(x, 2) for x in tf.shift_or_mirror_into_invertible_domain([1.52, -12.2, -0.2, 2, 4, 10.4])] 
1214      >>> assert res == [1.52, 9.2, 2.0, 2.0, 4.0, 10.4] 
1215      >>> tmp = tf([1])  # call with lower dimension 
1216   
1217      """ 
1218 -    def __init__(self, bounds): 
1219          """``x`` is defined in ``[lb - 3*al, ub + au + r - 2*al]`` with ``r = ub - lb + al + au``, 
1220          and ``x == transformation(x)`` in ``[lb + al, ub - au]``. 
1221          ``beta*x - alphal = beta*x - alphau`` is then defined in ``[lb, ub]``, 
1222   
1223          ``alphal`` and ``alphau`` represent the same value, but respectively numerically 
1224          better suited for values close to lb and ub. 
1225   
1226          """ 
1227          # super().__init__(bounds) # only in Python 3.x available 
1228          BoxConstraintsTransformationBase.__init__(self, bounds) 
1229          # super(BB, self).__init__(bounds) # is supposed to call initialize 
1230          # ## super(BoxConstraintsTransformationBase, self).__init__(bounds) # is probably invalid 
1231   
1232 -    def initialize(self, length=None): 
1233          """see ``__init__``""" 
1234          if length is None: 
1235              length = len(self.bounds) 
1236          max_i = min((len(self.bounds) - 1, length - 1)) 
1237          self._lb = array([self.bounds[min((i, max_i))][0] 
1238                            if self.bounds[min((i, max_i))][0] is not None else -np.Inf 
1239                            for i in xrange(length)], copy=False) 
1240          self._ub = array([self.bounds[min((i, max_i))][1] 
1241                            if self.bounds[min((i, max_i))][1] is not None else np.Inf 
1242                            for i in xrange(length)], copy=False) 
1243          lb = self._lb 
1244          ub = self._ub 
1245          # define added values for lower and upper bound 
1246          self._al = array([min([(ub[i] - lb[i]) / 2, (1 + np.abs(lb[i])) / 20]) 
1247                               if np.isfinite(lb[i]) else 1 for i in rglen(lb)], copy=False) 
1248          self._au = array([min([(ub[i] - lb[i]) / 2, (1 + np.abs(ub[i])) / 20]) 
1249                               if np.isfinite(ub[i]) else 1 for i in rglen(ub)], copy=False) 
1250   
1251 -    def __call__(self, solution_genotype, copy_if_changed=True, copy_always=False): 
1252          # about four times faster version of array([self._transform_i(x, i) for i, x in enumerate(solution_genotype)]) 
1253          # still, this makes a typical run on a test function two times slower, but there might be one too many copies 
1254          # during the transformations in gp 
1255          if len(self._lb) != len(solution_genotype): 
1256              self.initialize(len(solution_genotype)) 
1257          lb = self._lb 
1258          ub = self._ub 
1259          al = self._al 
1260          au = self._au 
1261   
1262          if copy_always or not isinstance(solution_genotype[0], float): 
1263              # transformed value is likely to be a float 
1264              y = np.array(solution_genotype, copy=True, dtype=float) 
1265              # if solution_genotype is not a float, copy value is disregarded 
1266              copy = False 
1267          else: 
1268              y = solution_genotype 
1269              copy = copy_if_changed 
1270          idx = (y < lb - 2 * al - (ub - lb) / 2.0) | (y > ub + 2 * au + (ub - lb) / 2.0) 
1271          if idx.any(): 
1272              r = 2 * (ub[idx] - lb[idx] + al[idx] + au[idx])  # period 
1273              s = lb[idx] - 2 * al[idx] - (ub[idx] - lb[idx]) / 2.0  # start 
1274              if copy: 
1275                  y = np.array(y, copy=True) 
1276                  copy = False 
1277              y[idx] -= r * ((y[idx] - s) // r)  # shift 
1278          idx = y > ub + au 
1279          if idx.any(): 
1280              if copy: 
1281                  y = np.array(y, copy=True) 
1282                  copy = False 
1283              y[idx] -= 2 * (y[idx] - ub[idx] - au[idx]) 
1284          idx = y < lb - al 
1285          if idx.any(): 
1286              if copy: 
1287                  y = np.array(y, copy=True) 
1288                  copy = False 
1289              y[idx] += 2 * (lb[idx] - al[idx] - y[idx]) 
1290          idx = y < lb + al 
1291          if idx.any(): 
1292              if copy: 
1293                  y = np.array(y, copy=True) 
1294                  copy = False 
1295              y[idx] = lb[idx] + (y[idx] - (lb[idx] - al[idx]))**2 / 4 / al[idx] 
1296          idx = y > ub - au 
1297          if idx.any(): 
1298              if copy: 
1299                  y = np.array(y, copy=True) 
1300                  copy = False 
1301              y[idx] = ub[idx] - (y[idx] - (ub[idx] + au[idx]))**2 / 4 / au[idx] 
1302          # assert Mh.vequals_approximately(y, BoxConstraintsTransformationBase.__call__(self, solution_genotype)) 
1303          return y 
1304      __call__.doc = BoxConstraintsTransformationBase.__doc__ 
1305      transform = __call__ 
1306 -    def idx_infeasible(self, solution_genotype): 
1307          """return indices of "infeasible" variables, that is, 
1308          variables that do not directly map into the feasible domain such that 
1309          ``tf.inverse(tf(x)) == x``. 
1310   
1311          """ 
1312          res = [i for i, x in enumerate(solution_genotype) if not self.is_feasible_i(x, i)] 
1313          return res 
1314 -    def is_feasible_i(self, x, i): 
1315          """return True if value ``x`` is in the invertible domain of 
1316          variable ``i`` 
1317   
1318          """ 
1319          lb = self._lb[self._index(i)] 
1320          ub = self._ub[self._index(i)] 
1321          al = self._al[self._index(i)] 
1322          au = self._au[self._index(i)] 
1323          return lb - al < x < ub + au 
1324 -    def is_loosely_feasible_i(self, x, i): 
1325          """never used""" 
1326          lb = self._lb[self._index(i)] 
1327          ub = self._ub[self._index(i)] 
1328          al = self._al[self._index(i)] 
1329          au = self._au[self._index(i)] 
1330          return lb - 2 * al - (ub - lb) / 2.0 <= x <= ub + 2 * au + (ub - lb) / 2.0 
1331   
1332 -    def shift_or_mirror_into_invertible_domain(self, solution_genotype, copy=False): 
1333          """Details: input ``solution_genotype`` is changed. The domain is 
1334          [lb - al, ub + au] and in [lb - 2*al - (ub - lb) / 2, lb - al] 
1335          mirroring is applied. 
1336   
1337          """ 
1338          assert solution_genotype is not None 
1339          if copy: 
1340              y = [val for val in solution_genotype] 
1341          else: 
1342              y = solution_genotype 
1343          if isinstance(y, np.ndarray) and not isinstance(y[0], float): 
1344              y = array(y, dtype=float) 
1345          for i in rglen(y): 
1346              lb = self._lb[self._index(i)] 
1347              ub = self._ub[self._index(i)] 
1348              al = self._al[self._index(i)] 
1349              au = self._au[self._index(i)] 
1350              # x is far from the boundary, compared to ub - lb 
1351              if y[i] < lb - 2 * al - (ub - lb) / 2.0 or y[i] > ub + 2 * au + (ub - lb) / 2.0: 
1352                  r = 2 * (ub - lb + al + au)  # period 
1353                  s = lb - 2 * al - (ub - lb) / 2.0  # start 
1354                  y[i] -= r * ((y[i] - s) // r)  # shift 
1355              if y[i] > ub + au: 
1356                  y[i] -= 2 * (y[i] - ub - au) 
1357              if y[i] < lb - al: 
1358                  y[i] += 2 * (lb - al - y[i]) 
1359          return y 
1360      shift_or_mirror_into_invertible_domain.__doc__ = BoxConstraintsTransformationBase.shift_or_mirror_into_invertible_domain.__doc__ + shift_or_mirror_into_invertible_domain.__doc__ 
1361   
1362 -    def _shift_or_mirror_into_invertible_i(self, x, i): 
1363          """shift into the invertible domain [lb - ab, ub + au], mirror close to 
1364          boundaries in order to get a smooth transformation everywhere 
1365   
1366          """ 
1367          assert x is not None 
1368          lb = self._lb[self._index(i)] 
1369          ub = self._ub[self._index(i)] 
1370          al = self._al[self._index(i)] 
1371          au = self._au[self._index(i)] 
1372          # x is far from the boundary, compared to ub - lb 
1373          if x < lb - 2 * al - (ub - lb) / 2.0 or x > ub + 2 * au + (ub - lb) / 2.0: 
1374              r = 2 * (ub - lb + al + au)  # period 
1375              s = lb - 2 * al - (ub - lb) / 2.0  # start 
1376              x -= r * ((x - s) // r)  # shift 
1377          if x > ub + au: 
1378              x -= 2 * (x - ub - au) 
1379          if x < lb - al: 
1380              x += 2 * (lb - al - x) 
1381          return x 
1382 -    def _transform_i(self, x, i): 
1383          """return transform of x in component i""" 
1384          x = self._shift_or_mirror_into_invertible_i(x, i) 
1385          lb = self._lb[self._index(i)] 
1386          ub = self._ub[self._index(i)] 
1387          al = self._al[self._index(i)] 
1388          au = self._au[self._index(i)] 
1389          if x < lb + al: 
1390              return lb + (x - (lb - al))**2 / 4 / al 
1391          elif x < ub - au: 
1392              return x 
1393          elif x < ub + 3 * au: 
1394              return ub - (x - (ub + au))**2 / 4 / au 
1395          else: 
1396              assert False  # shift removes this case 
1397              return ub + au - (x - (ub + au)) 
1398 -    def _inverse_i(self, y, i): 
1399          """return inverse of y in component i""" 
1400          lb = self._lb[self._index(i)] 
1401          ub = self._ub[self._index(i)] 
1402          al = self._al[self._index(i)] 
1403          au = self._au[self._index(i)] 
1404          if 1 < 3: 
1405              if not lb <= y <= ub: 
1406                  raise ValueError('argument of inverse must be within the given bounds') 
1407          if y < lb + al: 
1408              return (lb - al) + 2 * (al * (y - lb))**0.5 
1409          elif y < ub - au: 
1410              return y 
1411          else: 
1412              return (ub + au) - 2 * (au * (ub - y))**0.5 
1413   
1414  # ____________________________________________________________ 
1415  # ____________________________________________________________ 
1416  # 
1417 -class GenoPheno(object): 
1418      """Genotype-phenotype transformation. 
1419   
1420      Method `pheno` provides the transformation from geno- to phenotype, 
1421      that is from the internal representation to the representation used 
1422      in the objective function. Method `geno` provides the "inverse" pheno- 
1423      to genotype transformation. The geno-phenotype transformation comprises, 
1424      in this order: 
1425   
1426         - insert fixed variables (with the phenotypic and therefore quite 
1427           possibly "wrong" values) 
1428         - affine linear transformation (scaling and shift) 
1429         - user-defined transformation 
1430         - repair (e.g. into feasible domain due to boundaries) 
1431         - assign fixed variables their original phenotypic value 
1432   
1433      By default all transformations are the identity. The repair is only applied, 
1434      if the transformation is given as argument to the method `pheno`. 
1435   
1436      ``geno`` is only necessary, if solutions have been injected. 
1437   
1438      """ 
1439 -    def __init__(self, dim, scaling=None, typical_x=None, fixed_values=None, tf=None): 
1440          """return `GenoPheno` instance with fixed phenotypic dimension `dim`. 
1441   
1442          Keyword Arguments 
1443          ----------------- 
1444              `scaling` 
1445                  the diagonal of a scaling transformation matrix, multipliers 
1446                  in the genotyp-phenotyp transformation, see `typical_x` 
1447              `typical_x` 
1448                  ``pheno = scaling*geno + typical_x`` 
1449              `fixed_values` 
1450                  a dictionary of variable indices and values, like ``{0:2.0, 2:1.1}``, 
1451                  that are not subject to change, negative indices are ignored 
1452                  (they act like incommenting the index), values are phenotypic 
1453                  values. 
1454              `tf` 
1455                  list of two user-defined transformation functions, or `None`. 
1456   
1457                  ``tf[0]`` is a function that transforms the internal representation 
1458                  as used by the optimizer into a solution as used by the 
1459                  objective function. ``tf[1]`` does the back-transformation. 
1460                  For example :: 
1461   
1462                      tf_0 = lambda x: [xi**2 for xi in x] 
1463                      tf_1 = lambda x: [abs(xi)**0.5 fox xi in x] 
1464   
1465                  or "equivalently" without the `lambda` construct :: 
1466   
1467                      def tf_0(x): 
1468                          return [xi**2 for xi in x] 
1469                      def tf_1(x): 
1470                          return [abs(xi)**0.5 fox xi in x] 
1471   
1472                  ``tf=[tf_0, tf_1]`` is a reasonable way to guaranty that only positive 
1473                  values are used in the objective function. 
1474   
1475          Details 
1476          ------- 
1477          If ``tf_1`` is ommitted, the initial x-value must be given as genotype (because 
1478          the phenotype-genotype transformation is unknown in this case) and "injection" of 
1479          solutions might lead to unexpected results. 
1480   
1481          """ 
1482          self.N = dim 
1483          self.fixed_values = fixed_values 
1484          if tf is not None: 
1485              self.tf_pheno = tf[0] 
1486              self.tf_geno = tf[1]  # TODO: should not necessarily be needed 
1487              # r = np.random.randn(dim) 
1488              # assert all(tf[0](tf[1](r)) - r < 1e-7) 
1489              # r = np.random.randn(dim) 
1490              # assert all(tf[0](tf[1](r)) - r > -1e-7) 
1491              print("WARNING in class GenoPheno: user defined transformations have not been tested thoroughly") 
1492          else: 
1493              self.tf_geno = None 
1494              self.tf_pheno = None 
1495   
1496          if fixed_values: 
1497              if not isinstance(fixed_values, dict): 
1498                  raise _Error("fixed_values must be a dictionary {index:value,...}") 
1499              if max(fixed_values.keys()) >= dim: 
1500                  raise _Error("max(fixed_values.keys()) = " + str(max(fixed_values.keys())) + 
1501                      " >= dim=N=" + str(dim) + " is not a feasible index") 
1502              # convenience commenting functionality: drop negative keys 
1503              for k in list(fixed_values.keys()): 
1504                  if k < 0: 
1505                      fixed_values.pop(k) 
1506   
1507          def vec_is_default(vec, default_val=0): 
1508              """return True if `vec` has the value `default_val`, 
1509              None or [None] are also recognized as default 
1510   
1511              """ 
1512              # TODO: rather let default_val be a list of default values, cave comparison of arrays 
1513              try: 
1514                  if len(vec) == 1: 
1515                      vec = vec[0]  # [None] becomes None and is always default 
1516              except TypeError: 
1517                  pass  # vec is a scalar 
1518   
1519              if vec is None or vec == default_val: 
1520                  return True 
1521              try: 
1522                  if vec == array(None): 
1523                      return True 
1524              except NotImplementedError: 
1525                  pass 
1526              return False 
1527   
1528          self.scales = array(scaling) if scaling is not None else None 
1529          if vec_is_default(self.scales, 1): 
1530              self.scales = 1  # CAVE: 1 is not array(1) 
1531          elif self.scales.shape is not () and len(self.scales) != self.N: 
1532              raise _Error('len(scales) == ' + str(len(self.scales)) + 
1533                           ' does not match dimension N == ' + str(self.N)) 
1534   
1535          self.typical_x = array(typical_x) if typical_x is not None else None 
1536          if vec_is_default(self.typical_x, 0): 
1537              self.typical_x = 0 
1538          elif self.typical_x.shape is not () and len(self.typical_x) != self.N: 
1539              raise _Error('len(typical_x) == ' + str(len(self.typical_x)) + 
1540                           ' does not match dimension N == ' + str(self.N)) 
1541   
1542          if (self.scales is 1 and 
1543                  self.typical_x is 0 and 
1544                  self.fixed_values is None and 
1545                  self.tf_pheno is None): 
1546              self.isidentity = True 
1547          else: 
1548              self.isidentity = False 
1549   
1550 -    def pheno(self, x, into_bounds=(lambda x, copy=False: x if not copy else array(x, copy=copy)), 
1551                copy=True, copy_always=False, 
1552                archive=None, iteration=None): 
1553          """maps the genotypic input argument into the phenotypic space, see 
1554          help for class `GenoPheno` 
1555   
1556          Details 
1557          ------- 
1558          If ``copy``, values from ``x`` are copied if changed under the transformation. 
1559   
1560          """ 
1561          # TODO: copy_always seems superfluous, as it could be done in the calling code 
1562          if copy_always and not copy: 
1563              raise ValueError('arguments copy_always=' + str(copy_always) + 
1564                               ' and copy=' + str(copy) + ' have inconsistent values') 
1565          if copy_always: 
1566              x = array(x, copy=True) 
1567              copy = False 
1568   
1569          if self.isidentity: 
1570              y = into_bounds(x) # was into_bounds(x, False) before (bug before v0.96.22) 
1571          else: 
1572              if self.fixed_values is None: 
1573                  y = array(x, copy=copy)  # make a copy, in case 
1574              else:  # expand with fixed values 
1575                  y = list(x)  # is a copy 
1576                  for i in sorted(self.fixed_values.keys()): 
1577                      y.insert(i, self.fixed_values[i]) 
1578                  y = array(y, copy=False) 
1579              copy = False 
1580   
1581              if self.scales is not 1:  # just for efficiency 
1582                  y *= self.scales 
1583   
1584              if self.typical_x is not 0: 
1585                  y += self.typical_x 
1586   
1587              if self.tf_pheno is not None: 
1588                  y = array(self.tf_pheno(y), copy=False) 
1589   
1590              y = into_bounds(y, copy)  # copy is False 
1591   
1592              if self.fixed_values is not None: 
1593                  for i, k in list(self.fixed_values.items()): 
1594                      y[i] = k 
1595   
1596          if archive is not None: 
1597              archive.insert(y, geno=x, iteration=iteration) 
1598          return y 
1599   
1600 -    def geno(self, y, from_bounds=lambda x: x, 
1601               copy_if_changed=True, copy_always=False, 
1602               repair=None, archive=None): 
1603          """maps the phenotypic input argument into the genotypic space, 
1604          that is, computes essentially the inverse of ``pheno``. 
1605   
1606          By default a copy is made only to prevent to modify ``y``. 
1607   
1608          The inverse of the user-defined transformation (if any) 
1609          is only needed if external solutions are injected, it is not 
1610          applied to the initial solution x0. 
1611   
1612          Details 
1613          ======= 
1614          ``geno`` searches first in ``archive`` for the genotype of 
1615          ``y`` and returns the found value, typically unrepaired. 
1616          Otherwise, first ``from_bounds`` is applied, to revert a 
1617          projection into the bound domain, (if necessary) and ``pheno`` 
1618          is reverted. ``repair`` is applied last, and is usually the 
1619          method ``CMAEvolutionStrategy.repair_genotype`` that limits the 
1620          Mahalanobis norm of ``geno(y) - mean``. 
1621   
1622          """ 
1623          if archive is not None: 
1624              try: 
1625                  x = archive[y]['geno'] 
1626              except KeyError: 
1627                  x = None 
1628              if x is not None: 
1629                  if archive[y]['iteration'] < archive.last_iteration and repair is not None: 
1630                      x = repair(x, copy_if_changed=copy_always) 
1631                  return x 
1632          x = y 
1633          if copy_always: 
1634              x = array(y, copy=True) 
1635              copy = False 
1636          else: 
1637              copy = copy_if_changed 
1638   
1639          x = from_bounds(x)  # TODO should also take copy? 
1640   
1641          if self.isidentity: 
1642              if repair is not None: 
1643                  x = repair(x, copy) 
1644              return x 
1645   
1646          if copy:  # could be improved? 
1647              x = array(x, copy=True) 
1648              copy = False 
1649   
1650          # user-defined transformation 
1651          if self.tf_geno is not None: 
1652              x = array(self.tf_geno(x), copy=False) 
1653          elif self.tf_pheno is not None: 
1654              raise ValueError('t1 of options transformation was not defined but is needed as being the inverse of t0') 
1655   
1656          # affine-linear transformation: shift and scaling 
1657          if self.typical_x is not 0: 
1658              x -= self.typical_x 
1659          if self.scales is not 1:  # just for efficiency 
1660              x /= self.scales 
1661   
1662          # kick out fixed_values 
1663          if self.fixed_values is not None: 
1664              # keeping the transformed values does not help much 
1665              # therefore it is omitted 
1666              if 1 < 3: 
1667                  keys = sorted(self.fixed_values.keys()) 
1668                  x = array([x[i] for i in range(len(x)) if i not in keys], copy=False) 
1669          # repair injected solutions 
1670          if repair is not None: 
1671              x = repair(x, copy) 
1672          return x 
1673   
1674  # ____________________________________________________________ 
1675  # ____________________________________________________________ 
1676  # check out built-in package abc: class ABCMeta, abstractmethod, abstractproperty... 
1677  # see http://docs.python.org/whatsnew/2.6.html PEP 3119 abstract base classes 
1678  # 
1679 -class OOOptimizer(object): 
1680      """"abstract" base class for an OO optimizer interface. 
1681   
1682       Relevant methods are `__init__`, `ask`, `tell`, `stop`, `result`, 
1683       and `optimize`. Only `optimize` is fully implemented in this base 
1684       class. 
1685   
1686      Examples 
1687      -------- 
1688      All examples minimize the function `elli`, the output is not shown. 
1689      (A preferred environment to execute all examples is ``ipython -pylab``.) 
1690      First we need :: 
1691   
1692          from cma import CMAEvolutionStrategy, CMADataLogger  # CMAEvolutionStrategy derives from the OOOptimizer class 
1693          elli = lambda x: sum(1e3**((i-1.)/(len(x)-1.)*x[i])**2 for i in range(len(x))) 
1694   
1695      The shortest example uses the inherited method `OOOptimizer.optimize()`:: 
1696   
1697          res = CMAEvolutionStrategy(8 * [0.1], 0.5).optimize(elli) 
1698   
1699      The input parameters to `CMAEvolutionStrategy` are specific to this 
1700      inherited class. The remaining functionality is based on interface 
1701      defined by `OOOptimizer`. We might have a look at the result:: 
1702   
1703          print(res[0])  # best solution and 
1704          print(res[1])  # its function value 
1705   
1706      `res` is the return value from method `CMAEvolutionStrategy.result()` 
1707      appended with CMAEvolutionStrategy.logger. 
1708      In order to display more exciting output we do :: 
1709   
1710          res[-1].plot()  # if matplotlib is available 
1711   
1712      Virtually the same example can be written with an explicit loop 
1713      instead of using `optimize()`. This gives the necessary insight into 
1714      the `OOOptimizer` class interface and gives entire control over the 
1715      iteration loop:: 
1716   
1717          optim = CMAEvolutionStrategy(9 * [0.5], 0.3)  # a new CMAEvolutionStrategy instance calling CMAEvolutionStrategy.__init__() 
1718          logger = CMADataLogger(optim)  # get a logger instance, we could also use the instance optim.logger 
1719   
1720          # this loop resembles optimize() 
1721          while not optim.stop(): # iterate 
1722              X = optim.ask()     # get candidate solutions 
1723              f = [elli(x) for x in X]  # evaluate solutions 
1724              #  maybe do something else that needs to be done 
1725              optim.tell(X, f)    # do all the real work: prepare for next iteration 
1726              optim.disp(20)      # display info every 20th iteration 
1727              logger.add()        # log another "data line" 
1728   
1729          # final output 
1730          print('termination by', optim.stop()) 
1731          print('best f-value =', optim.result()[1]) 
1732          print('best solution =', optim.result()[0]) 
1733          logger.plot()  # if matplotlib is available 
1734          raw_input('press enter to continue')  # prevents exiting and closing figures 
1735   
1736      Details 
1737      ------- 
1738      Most of the work is done in the method `tell(...)`. The method `result()` returns 
1739      more useful output. 
1740   
1741      """ 
1742 -    def __init__(self, xstart, **more_args): 
1743          """``xstart`` is a mandatory argument""" 
1744          self.xstart = xstart 
1745          self.more_args = more_args 
1746          self.initialize() 
1747 -    def initialize(self): 
1748          """(re-)set to the initial state""" 
1749          self.countiter = 0 
1750          self.xcurrent = self.xstart[:] 
1751          raise NotImplementedError('method initialize() must be implemented in derived class') 
1752 -    def ask(self): 
1753          """abstract method, AKA "get" or "sample_distribution", deliver 
1754          new candidate solution(s), a list of "vectors" 
1755          """ 
1756          raise NotImplementedError('method ask() must be implemented in derived class') 
1757 -    def tell(self, solutions, function_values): 
1758          """abstract method, AKA "update", prepare for next iteration""" 
1759          self.countiter += 1 
1760          raise NotImplementedError('method tell() must be implemented in derived class') 
1761 -    def stop(self): 
1762          """abstract method, return satisfied termination conditions in 
1763          a dictionary like ``{'termination reason': value, ...}``, 
1764          for example ``{'tolfun': 1e-12}``, or the empty dictionary ``{}``. 
1765          The implementation of `stop()` should prevent an infinite 
1766          loop. 
1767          """ 
1768          raise NotImplementedError('method stop() is not implemented') 
1769 -    def disp(self, modulo=None): 
1770          """abstract method, display some iteration infos if 
1771          ``self.iteration_counter % modulo == 0`` """ 
1772          pass  # raise NotImplementedError('method disp() is not implemented') 
1773 -    def result(self): 
1774          """abstract method, return ``(x, f(x), ...)``, that is, the 
1775          minimizer, its function value, ...""" 
1776          raise NotImplementedError('method result() is not implemented') 
1777   
1778      # previous ordering: 
1779      #    def optimize(self, objectivefct, 
1780      #                 logger=None, verb_disp=20, 
1781      #                 iterations=None, min_iterations=1, 
1782      #                 call_back=None): 
1783 -    def optimize(self, objective_fct, 
1784                   iterations=None, min_iterations=1, args=(), 
1785                   verb_disp=None, logger=None, 
1786                   call_back=None): 
1787          """find minimizer of `objective_fct`. 
1788   
1789          CAVEAT: the return value for `optimize` is versatile and might 
1790          change in near future to ``self``. 
1791   
1792          Arguments 
1793          --------- 
1794   
1795              `objective_fct` 
1796                  function be to minimized 
1797              `iterations` 
1798                  number of (maximal) iterations, while ``not self.stop()`` 
1799              `min_iterations` 
1800                  minimal number of iterations, even if ``not self.stop()`` 
1801              `args` 
1802                  arguments passed to `objective_fct` 
1803              `verb_disp` 
1804                  print to screen every `verb_disp` iteration, if ``None`` 
1805                  the value from ``self.logger`` is "inherited", if 
1806                  available. 
1807              ``logger`` 
1808                  a `BaseDataLogger` instance, which must be compatible 
1809                  with the type of ``self``. 
1810              ``call_back`` 
1811                  call back function called like ``call_back(self)`` or 
1812                  a list of call back functions. 
1813   
1814          ``return self.result() + (self.stop(), self, logger)`` which 
1815          might change in near future. 
1816   
1817          Example 
1818          ------- 
1819          >>> import cma 
1820          >>> res = cma.CMAEvolutionStrategy(7 * [0.1], 0.5).optimize(cma.fcts.rosen, verb_disp=100) 
1821          (4_w,9)-CMA-ES (mu_w=2.8,w_1=49%) in dimension 7 (seed=630721393) 
1822          Iterat #Fevals   function value     axis ratio  sigma   minstd maxstd min:sec 
1823              1       9 3.163954777181882e+01 1.0e+00 4.12e-01  4e-01  4e-01 0:0.0 
1824              2      18 3.299006223906629e+01 1.0e+00 3.60e-01  3e-01  4e-01 0:0.0 
1825              3      27 1.389129389866704e+01 1.1e+00 3.18e-01  3e-01  3e-01 0:0.0 
1826            100     900 2.494847340045985e+00 8.6e+00 5.03e-02  2e-02  5e-02 0:0.3 
1827            200    1800 3.428234862999135e-01 1.7e+01 3.77e-02  6e-03  3e-02 0:0.5 
1828            300    2700 3.216640032470860e-04 5.6e+01 6.62e-03  4e-04  9e-03 0:0.8 
1829            400    3600 6.155215286199821e-12 6.6e+01 7.44e-06  1e-07  4e-06 0:1.1 
1830            438    3942 1.187372505161762e-14 6.0e+01 3.27e-07  4e-09  9e-08 0:1.2 
1831            438    3942 1.187372505161762e-14 6.0e+01 3.27e-07  4e-09  9e-08 0:1.2 
1832          ('termination by', {'tolfun': 1e-11}) 
1833          ('best f-value =', 1.1189867885201275e-14) 
1834          ('solution =', array([ 1.        ,  1.        ,  1.        ,  0.99999999,  0.99999998, 
1835                  0.99999996,  0.99999992])) 
1836          >>> print(res[0]) 
1837          [ 1.          1.          1.          0.99999999  0.99999998  0.99999996 
1838            0.99999992] 
1839   
1840          """ 
1841          assert iterations is None or min_iterations <= iterations 
1842          if not hasattr(self, 'logger'): 
1843              self.logger = logger 
1844          if logger is None: 
1845              logger = self.logger 
1846          self.logger_from_optimize_method_call = logger 
1847          if not isinstance(call_back, list): 
1848              call_back = [call_back] 
1849   
1850          citer = 0 
1851          while not self.stop() or citer < min_iterations: 
1852              if iterations is not None and citer >= iterations: 
1853                  return self.result() 
1854              citer += 1 
1855   
1856              X = self.ask()  # deliver candidate solutions 
1857              fitvals = [objective_fct(x, *args) for x in X] 
1858              self.tell(X, fitvals)  # all the work is done here 
1859              self.disp(verb_disp) 
1860              for f in call_back: 
1861                  if f is not None: 
1862                      f(self) 
1863              logger.add(self) if logger else None 
1864   
1865          # signal logger that we left the loop 
1866          # TODO: this is very ugly, because it assumes modulo keyword 
1867          #       argument *and* modulo attribute to be available 
1868          try: 
1869              logger.add(self, modulo=bool(logger.modulo)) if logger else None 
1870          except TypeError: 
1871              print('  suppressing the final call of the logger in OOOptimizer.optimize (modulo keyword parameter not available)') 
1872          except AttributeError: 
1873              print('  suppressing the final call of the logger in OOOptimizer.optimize (modulo attribute not available)') 
1874          if verb_disp: 
1875              self.disp(1) 
1876          if verb_disp in (1, True): 
1877              print('termination by', self.stop()) 
1878              print('best f-value =', self.result()[1]) 
1879              print('solution =', self.result()[0]) 
1880   
1881          return self.result() + (self.stop(), self, logger) 
1882   
1883 -class CMAAdaptSigmaBase(object): 
1884      """step-size adaptation base class, implementing hsig functionality 
1885      via an isotropic evolution path. 
1886   
1887      """ 
1888 -    def __init__(self, *args, **kwargs): 
1889          self.is_initialized_base = False 
1890          self._ps_updated_iteration = -1 
1891 -    def initialize_base(self, es): 
1892          """set parameters and state variable based on dimension, 
1893          mueff and possibly further options. 
1894   
1895          """ 
1896          self.cs = (es.sp.mueff + 2) / (es.N + es.sp.mueff + 3) 
1897          self.ps = np.zeros(es.N) 
1898          self.is_initialized_base = True 
1899 -    def _update_ps(self, es): 
1900          """update the isotropic evolution path 
1901   
1902          :type es: CMAEvolutionStrategy 
1903          """ 
1904          if not self.is_initialized_base: 
1905              self.initialize_base(es) 
1906          if self._ps_updated_iteration == es.countiter: 
1907              return 
1908          if es.countiter <= es.itereigenupdated: 
1909              # es.B and es.D must/should be those from the last iteration 
1910              assert es.countiter >= es.itereigenupdated 
1911              _print_warning('distribution transformation (B and D) have been updated before ps could be computed', 
1912                            '_update_ps', 'CMAAdaptSigmaBase') 
1913          z = dot(es.B, (1. / es.D) * dot(es.B.T, (es.mean - es.mean_old) / es.sigma_vec)) 
1914          z *= es.sp.mueff**0.5 / es.sigma / es.sp.cmean 
1915          # self.cs or es.sp.cs could be used here 
1916          self.ps = (1 - self.cs) * self.ps + sqrt(self.cs * (2 - self.cs)) * z 
1917          self._ps_updated_iteration = es.countiter 
1918 -    def hsig(self, es): 
1919          """return "OK-signal" for rank-one update, `True` (OK) or `False` 
1920          (stall rank-one update), based on the length of an evolution path 
1921   
1922          """ 
1923          self._update_ps(es) 
1924          if self.ps is None: 
1925              return True 
1926          squared_sum = sum(self.ps**2) / (1 - (1 - self.cs)**(2 * es.countiter)) 
1927          # correction with self.countiter seems not necessary, 
1928          # as pc also starts with zero 
1929          return squared_sum / es.N - 1 < 1 + 4. / (es.N + 1) 
1930 -    def update(self, es, **kwargs): 
1931          """update ``es.sigma``""" 
1932          self._update_ps(es) 
1933          raise NotImplementedError('must be implemented in a derived class') 
1934 -class CMAAdaptSigmaNone(CMAAdaptSigmaBase): 
1935 -    def update(self, es, **kwargs): 
1936          """no update, ``es.sigma`` remains constant. 
1937   
1938          :param es: ``CMAEvolutionStrategy`` class instance 
1939          :param kwargs: whatever else is needed to update ``es.sigma`` 
1940   
1941          """ 
1942          pass 
1943 -class CMAAdaptSigmaDistanceProportional(CMAAdaptSigmaBase): 
1944      """artificial setting of ``sigma`` for test purposes, e.g. 
1945      to simulate optimal progress rates. 
1946   
1947      """ 
1948 -    def __init__(self, coefficient=1.2): 
1949          self.coefficient = coefficient 
1950 -    def update(self, es, **kwargs): 
1951          # optimal step-size is 
1952          es.sigma = self.coefficient * self.sp.mueff * sum(self.mean**2)**0.5 / self.N 
1953 -class CMAAdaptSigmaCSA(CMAAdaptSigmaBase): 
1954 -    def __init__(self): 
1955          """postpone initialization to a method call where dimension and mueff should be known. 
1956   
1957          """ 
1958          self.is_initialized = False 
1959 -    def initialize(self, es): 
1960          """set parameters and state variable based on dimension, 
1961          mueff and possibly further options. 
1962   
1963          """ 
1964          self.disregard_length_setting = True if es.opts['CSA_disregard_length'] else False 
1965          if es.opts['CSA_clip_length_value'] is not None: 
1966              try: 
1967                  if len(es.opts['CSA_clip_length_value']) == 0: 
1968                      es.opts['CSA_clip_length_value'] = [-np.Inf, np.Inf] 
1969                  elif len(es.opts['CSA_clip_length_value']) == 1: 
1970                      es.opts['CSA_clip_length_value'] = [-np.Inf, es.opts['CSA_clip_length_value'][0]] 
1971                  elif len(es.opts['CSA_clip_length_value']) == 2: 
1972                      es.opts['CSA_clip_length_value'] = np.sort(es.opts['CSA_clip_length_value']) 
1973                  else: 
1974                      raise ValueError('option CSA_clip_length_value should be a number of len(.) in [1,2]') 
1975              except TypeError:  # len(...) failed 
1976                  es.opts['CSA_clip_length_value'] = [-np.Inf, es.opts['CSA_clip_length_value']] 
1977              es.opts['CSA_clip_length_value'] = list(np.sort(es.opts['CSA_clip_length_value'])) 
1978              if es.opts['CSA_clip_length_value'][0] > 0 or es.opts['CSA_clip_length_value'][1] < 0: 
1979                  raise ValueError('option CSA_clip_length_value must be a single positive or a negative and a positive number') 
1980          self.cs = (es.sp.mueff + 2) / (es.N + es.sp.mueff + 3) 
1981          self.damps = es.opts['CSA_dampfac'] * (0.5 + 
1982                                            0.5 * min([1, (es.sp.lam_mirr / (0.159 * es.sp.popsize) - 1)**2])**1 + 
1983                                            2 * max([0, ((es.sp.mueff - 1) / (es.N + 1))**es.opts['CSA_damp_mueff_exponent'] - 1]) + 
1984                                            self.cs 
1985                                            ) 
1986          self.max_delta_log_sigma = 1  # in symmetric use (strict lower bound is -cs/damps anyway) 
1987   
1988          if self.disregard_length_setting: 
1989              es.opts['CSA_clip_length_value'] = [0, 0] 
1990              self.cs = (es.sp.mueff + 1)**0.5 / (es.N**0.5 + 2 * es.sp.mueff**0.5) 
1991              self.damps = es.opts['CSA_dampfac'] * 1  # * (1.1 - 1/(es.N+1)**0.5) 
1992          if es.opts['verbose'] > 1 or self.disregard_length_setting or 11 < 3: 
1993              print('SigmaCSA Parameters') 
1994              for k, v in self.__dict__.items(): 
1995                  print('  ', k, ':', v) 
1996          self.ps = np.zeros(es.N) 
1997          self._ps_updated_iteration = -1 
1998          self.is_initialized = True 
1999   
2000 -    def _update_ps(self, es): 
2001          if not self.is_initialized: 
2002              self.initialize(es) 
2003          if self._ps_updated_iteration == es.countiter: 
2004              return 
2005          z = dot(es.B, (1. / es.D) * dot(es.B.T, (es.mean - es.mean_old) / es.sigma_vec)) 
2006          z *= es.sp.mueff**0.5 / es.sigma / es.sp.cmean 
2007          # zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz 
2008          if es.opts['CSA_clip_length_value'] is not None: 
2009              vals = es.opts['CSA_clip_length_value'] 
2010              min_len = es.N**0.5 + vals[0] * es.N / (es.N + 2) 
2011              max_len = es.N**0.5 + vals[1] * es.N / (es.N + 2) 
2012              act_len = sum(z**2)**0.5 
2013              new_len = Mh.minmax(act_len, min_len, max_len) 
2014              if new_len != act_len: 
2015                  z *= new_len / act_len 
2016                  # z *= (es.N / sum(z**2))**0.5  # ==> sum(z**2) == es.N 
2017                  # z *= es.const.chiN / sum(z**2)**0.5 
2018          self.ps = (1 - self.cs) * self.ps + sqrt(self.cs * (2 - self.cs)) * z 
2019          self._ps_updated_iteration = es.countiter 
2020 -    def update(self, es, **kwargs): 
2021          self._update_ps(es)  # caveat: if es.B or es.D are already updated and ps is not, this goes wrong! 
2022          if es.opts['CSA_squared']: 
2023              s = (sum(self.ps**2) / es.N - 1) / 2 
2024              # sum(self.ps**2) / es.N has mean 1 and std sqrt(2/N) and is skewed 
2025              # divided by 2 to have the derivative d/dx (x**2 / N - 1) for x**2=N equal to 1 
2026          else: 
2027              s = sum(self.ps**2)**0.5 / es.const.chiN - 1 
2028          s *= self.cs / self.damps 
2029          s_clipped = Mh.minmax(s, -self.max_delta_log_sigma, self.max_delta_log_sigma) 
2030          es.sigma *= np.exp(s_clipped) 
2031          # "error" handling 
2032          if s_clipped != s: 
2033              _print_warning('sigma change exp(' + str(s) + ') = ' + str(np.exp(s)) + 
2034                            ' clipped to exp(+-' + str(self.max_delta_log_sigma) + ')', 
2035                            'update', 
2036                            'CMAAdaptSigmaCSA', 
2037                            es.countiter, es.opts['verbose']) 
2038 -class CMAAdaptSigmaMedianImprovement(CMAAdaptSigmaBase): 
2039      """Compares median fitness against a fitness percentile of the previous iteration, 
2040      see Ait ElHara et al, GECCO 2013. 
2041   
2042      """ 
2043 -    def __init__(self): 
2044          CMAAdaptSigmaBase.__init__(self) 
2045 -    def initialize(self, es): 
2046          r = es.sp.mueff / es.popsize 
2047          self.index_to_compare = 0.5 * (r**0.5 + 2.0 * (1 - r**0.5) / log(es.N + 9)**2) * (es.popsize)  # TODO 
2048          self.index_to_compare = (0.30 if not es.opts['vv'] else es.opts['vv']) * es.popsize  # TODO 
2049          self.damp = 2 - 2 / es.N  # sign-rule: 2 
2050          self.c = 0.3  # sign-rule needs <= 0.3 
2051          self.s = 0  # averaged statistics, usually between -1 and +1 
2052 -    def update(self, es, **kwargs): 
2053          if es.countiter < 2: 
2054              self.initialize(es) 
2055              self.fit = es.fit.fit 
2056          else: 
2057              ft1, ft2 = self.fit[int(self.index_to_compare)], self.fit[int(np.ceil(self.index_to_compare))] 
2058              ftt1, ftt2 = es.fit.fit[(es.popsize - 1) // 2], es.fit.fit[int(np.ceil((es.popsize - 1) / 2))] 
2059              pt2 = self.index_to_compare - int(self.index_to_compare) 
2060              # ptt2 = (es.popsize - 1) / 2 - (es.popsize - 1) // 2  # not in use 
2061              s = 0 
2062              if 1 < 3: 
2063                  s += (1 - pt2) * sum(es.fit.fit <= self.fit[int(np.ceil(self.index_to_compare))]) 
2064                  s += pt2 * sum(es.fit.fit < self.fit[int(self.index_to_compare)]) 
2065                  s -= (es.popsize + 1) / 2 
2066                  s *= 2 / es.popsize  # the range was popsize, is 2 
2067              self.s = (1 - self.c) * self.s + self.c * s 
2068              es.sigma *= exp(self.s / self.damp) 
2069          # es.more_to_write.append(10**(self.s)) 
2070   
2071          #es.more_to_write.append(10**((2 / es.popsize) * (sum(es.fit.fit < self.fit[int(self.index_to_compare)]) - (es.popsize + 1) / 2))) 
2072          # # es.more_to_write.append(10**(self.index_to_compare - sum(self.fit <= es.fit.fit[es.popsize // 2]))) 
2073          # # es.more_to_write.append(10**(np.sign(self.fit[int(self.index_to_compare)] - es.fit.fit[es.popsize // 2]))) 
2074          self.fit = es.fit.fit 
2075 -class CMAAdaptSigmaTPA(CMAAdaptSigmaBase): 
2076      """two point adaptation for step-size sigma. Relies on a specific 
2077      sampling of the first two offspring, whose objective function 
2078      value ranks are used to decide on the step-size change. 
2079   
2080      Example 
2081      ======= 
2082   
2083      >>> import cma 
2084      >>> cma.CMAOptions('adapt').pprint() 
2085      >>> es = cma.CMAEvolutionStrategy(10 * [0.2], 0.1, {'AdaptSigma': cma.CMAAdaptSigmaTPA, 'ftarget': 1e-8}) 
2086      >>> es.optimize(cma.fcts.rosen) 
2087      >>> assert 'ftarget' in es.stop() 
2088      >>> assert es.result()[1] <= 1e-8 
2089      >>> assert es.result()[2] < 6500  # typically < 5500 
2090   
2091      References: loosely based on Hansen 2008, CMA-ES with Two-Point 
2092      Step-Size Adaptation, more tightly based on an upcoming paper by 
2093      Hansen et al. 
2094   
2095      """ 
2096 -    def __init__(self, dimension=None, opts=None): 
2097          CMAAdaptSigmaBase.__init__(self) # base class provides method hsig() 
2098          self.initialized = False 
2099 -    def initialize(self, N): 
2100          self.sp = _BlancClass() 
2101          self.sp.damp = eval('N**0.5')  # why do we need 10 <-> exp(1/10) == 1.1? 2 should be fine!? 
2102          self.sp.dampup = 1.0 * self.sp.damp  # 0.5 fails to converge on the Rastrigin function 
2103          self.sp.dampdown = 1.0 * self.sp.damp 
2104          self.sp.c = 1.0  # rank difference is asymetric and therefore the switch from increase to decrease takes too long 
2105          self.sp.z_exponent = 0.5  # sign(z) * abs(z)**z_exponent, 0.5 seems better with larger popsize 
2106          self.sp.sigma_fac = 1.0  # (obsolete) 0.5 feels better, but no evidence whether it is 
2107          self.sp.relative_to_delta_mean = True  # (obsolete) 
2108          self.s = 0  # the state variable 
2109          self.last = None 
2110          self.initialized = True 
2111          return self 
2112 -    def update(self, es, function_values, **kwargs): 
2113          """the first and second value in ``function_values`` 
2114          must reflect two mirrored solutions sampled 
2115          in direction / in opposite direction of 
2116          the previous mean shift, respectively. 
2117   
2118          """ 
2119          # TODO: on the linear function, the two mirrored samples lead 
2120          # to a sharp increase of condition of the covariance matrix. 
2121          # They should not be used to update the covariance matrix, 
2122          # if the step-size inreases quickly. 
2123          if not self.initialized: 
2124              self.initialize(es.N) 
2125          if 1 < 3: 
2126              # use the ranking difference of the mirrors for adaptation 
2127              # damp = 5 should be fine 
2128              z = np.where(es.fit.idx == 1)[0][0] - np.where(es.fit.idx == 0)[0][0] 
2129              z /= es.popsize - 1 
2130          self.s = (1 - self.sp.c) * self.s + self.sp.c * np.sign(z) * np.abs(z)**self.sp.z_exponent 
2131          if self.s > 0: 
2132              es.sigma *= exp(self.s / self.sp.dampup) 
2133          else: 
2134              es.sigma *= exp(self.s / self.sp.dampdown) 
2135          #es.more_to_write.append(10**z) 
2136   
2137   
2138  # ____________________________________________________________ 
2139  # ____________________________________________________________ 
2140  # 
2141 -class CMAEvolutionStrategy(OOOptimizer): 
2142      """CMA-ES stochastic optimizer class with ask-and-tell interface. 
2143   
2144      Calling Sequences 
2145      ================= 
2146   
2147          ``es = CMAEvolutionStrategy(x0, sigma0)`` 
2148   
2149          ``es = CMAEvolutionStrategy(x0, sigma0, opts)`` 
2150   
2151          ``res = CMAEvolutionStrategy(x0, sigma0).optimize(objective_fct)`` 
2152   
2153      CAVEAT: return value of `optimize` might become ``optim`` in near 
2154      future. 
2155   
2156      Arguments 
2157      ========= 
2158          `x0` 
2159              initial solution, starting point. `x0` is given as "genotype" 
2160              which means, if:: 
2161   
2162                  opts={'transformation':[transform, inverse]} 
2163   
2164              is given (``inverse`` can be ``None``), then ``transform(x0)`` 
2165              is the "phenotypic" initial solution and 
2166              ``objective_function(transform(x0))`` is the objective 
2167              function value of ``x0``. 
2168   
2169          `sigma0` 
2170              initial standard deviation.  The problem variables should 
2171              have been scaled, such that a single standard deviation 
2172              on all variables is useful and the optimum is expected to 
2173              lie within about `x0` +- ``3*sigma0``. See also options 
2174              `scaling_of_variables`. Often one wants to check for 
2175              solutions close to the initial point. This allows, 
2176              for example, for an easier check of consistency of the 
2177              objective function and its interfacing with the optimizer. 
2178              In this case, a much smaller `sigma0` is advisable. 
2179          `opts` 
2180              options, a dictionary with optional settings, 
2181              see class `CMAOptions`. 
2182   
2183      Main interface / usage 
2184      ====================== 
2185      The interface is inherited from the generic `OOOptimizer` 
2186      class (see also there). An object instance is generated from 
2187   
2188          es = cma.CMAEvolutionStrategy(8 * [0.5], 0.2) 
2189   
2190      The least verbose interface is via the optimize method:: 
2191   
2192          es.optimize(objective_func) 
2193          res = es.result() 
2194   
2195      More verbosely, the optimization is done using the 
2196      methods ``stop``, ``ask``, and ``tell`` :: 
2197   
2198          while not es.stop(): 
2199              solutions = es.ask() 
2200              es.tell(solutions, [cma.fcts.rosen(s) for s in solutions]) 
2201   
2202   
2203      where ``ask`` delivers new candidate solutions and ``tell`` updates 
2204      the ``optim`` instance by passing the respective function values 
2205      (the objective function ``cma.fcts.rosen`` can be replaced by any 
2206      properly defined objective function, see ``cma.fcts`` for more 
2207      examples). 
2208   
2209      The class `CMAEvolutionStrategy` also provides:: 
2210   
2211          (solutions, func_values) = es.ask_and_eval(objective_func) 
2212   
2213      and an entire optimization can also be written like:: 
2214   
2215          while not es.stop(): 
2216              es.tell(*es.ask_and_eval(objective_func)) 
2217   
2218      Besides for termination criteria, in CMA-ES only the ranks of the 
2219      `func_values` are relevant. 
2220   
2221      Attributes and Properties 
2222      ========================= 
2223          - `inputargs` -- passed input arguments 
2224          - `inopts` -- passed options 
2225          - `opts` -- actually used options, some of them can be changed any 
2226            time, see class `CMAOptions` 
2227          - `popsize` -- population size lambda, number of candidate solutions 
2228            returned by `ask()` 
2229          - `logger` -- a `CMADataLogger` instance utilized by `optimize` 
2230   
2231      Examples 
2232      ======== 
2233      Super-short example, with output shown: 
2234   
2235      >>> import cma 
2236      >>> # construct an object instance in 4-D, sigma0=1: 
2237      >>> es = cma.CMAEvolutionStrategy(4 * [1], 1, {'seed':234}) 
2238      (4_w,8)-CMA-ES (mu_w=2.6,w_1=52%) in dimension 4 (seed=234) 
2239      >>> 
2240      >>> # optimize the ellipsoid function 
2241      >>> es.optimize(cma.fcts.elli, verb_disp=1) 
2242      Iterat #Fevals   function value     axis ratio  sigma   minstd maxstd min:sec 
2243          1       8 2.093015112685775e+04 1.0e+00 9.27e-01  9e-01  9e-01 0:0.0 
2244          2      16 4.964814235917688e+04 1.1e+00 9.54e-01  9e-01  1e+00 0:0.0 
2245          3      24 2.876682459926845e+05 1.2e+00 1.02e+00  9e-01  1e+00 0:0.0 
2246        100     800 6.809045875281943e-01 1.3e+02 1.41e-02  1e-04  1e-02 0:0.2 
2247        200    1600 2.473662150861846e-10 8.0e+02 3.08e-05  1e-08  8e-06 0:0.5 
2248        233    1864 2.766344961865341e-14 8.6e+02 7.99e-07  8e-11  7e-08 0:0.6 
2249      >>> 
2250      >>> cma.pprint(es.result()) 
2251      (array([ -1.98546755e-09,  -1.10214235e-09,   6.43822409e-11, 
2252              -1.68621326e-11]), 
2253       4.5119610261406537e-16, 
2254       1666, 
2255       1672, 
2256       209, 
2257       array([ -9.13545269e-09,  -1.45520541e-09,  -6.47755631e-11, 
2258              -1.00643523e-11]), 
2259       array([  3.20258681e-08,   3.15614974e-09,   2.75282215e-10, 
2260               3.27482983e-11])) 
2261      >>> assert es.result()[1] < 1e-9 
2262      >>> help(es.result) 
2263      Help on method result in module cma: 
2264   
2265      result(self) method of cma.CMAEvolutionStrategy instance 
2266          return ``(xbest, f(xbest), evaluations_xbest, evaluations, iterations, pheno(xmean), effective_stds)`` 
2267   
2268   
2269      The optimization loop can also be written explicitly. 
2270   
2271      >>> import cma 
2272      >>> es = cma.CMAEvolutionStrategy(4 * [1], 1) 
2273      >>> while not es.stop(): 
2274      ...    X = es.ask() 
2275      ...    es.tell(X, [cma.fcts.elli(x) for x in X]) 
2276      ...    es.disp() 
2277      <output omitted> 
2278   
2279      achieving the same result as above. 
2280   
2281      An example with lower bounds (at zero) and handling infeasible 
2282      solutions: 
2283   
2284      >>> import cma 
2285      >>> import numpy as np 
2286      >>> es = cma.CMAEvolutionStrategy(10 * [0.2], 0.5, {'bounds': [0, np.inf]}) 
2287      >>> while not es.stop(): 
2288      ...     fit, X = [], [] 
2289      ...     while len(X) < es.popsize: 
2290      ...         curr_fit = None 
2291      ...         while curr_fit in (None, np.NaN): 
2292      ...             x = es.ask(1)[0] 
2293      ...             curr_fit = cma.fcts.somenan(x, cma.fcts.elli) # might return np.NaN 
2294      ...         X.append(x) 
2295      ...         fit.append(curr_fit) 
2296      ...     es.tell(X, fit) 
2297      ...     es.logger.add() 
2298      ...     es.disp() 
2299      <output omitted> 
2300      >>> 
2301      >>> assert es.result()[1] < 1e-9 
2302      >>> assert es.result()[2] < 9000  # by internal termination 
2303      >>> # es.logger.plot()  # will plot data 
2304      >>> # cma.show()  # display plot window 
2305   
2306      An example with user-defined transformation, in this case to realize 
2307      a lowwer bound of 2. 
2308   
2309      >>> es = cma.CMAEvolutionStrategy(5 * [3], 1, 
2310      ...                 {"transformation": [lambda x: x**2+2, None]}) 
2311      >>> es.optimize(cma.fcts.rosen) 
2312      <output omitted> 
2313      >>> assert cma.fcts.rosen(es.result()[0]) < 1e-6 + 5.530760944396627e+02 
2314      >>> assert es.result()[2] < 3300 
2315   
2316      The inverse transformation is (only) necessary if the `BoundPenalty` 
2317      boundary handler is used at the same time. 
2318   
2319      The ``CMAEvolutionStrategy`` class also provides a default logger 
2320      (cave: files are overwritten when the logger is used with the same 
2321      filename prefix): 
2322   
2323      >>> import cma 
2324      >>> es = cma.CMAEvolutionStrategy(4 * [0.2], 0.5, {'verb_disp': 0}) 
2325      >>> es.logger.disp_header()  # to understand the print of disp 
2326      Iterat Nfevals  function value    axis ratio maxstd   minstd 
2327      >>> while not es.stop(): 
2328      ...     X = es.ask() 
2329      ...     es.tell(X, [cma.fcts.sphere(x) for x in X]) 
2330      ...     es.logger.add()  # log current iteration 
2331      ...     es.logger.disp([-1])  # display info for last iteration 
2332      1      8 2.72769793021748e+03 1.0e+00 4.05e-01 3.99e-01 
2333      2     16 6.58755537926063e+03 1.1e+00 4.00e-01 3.39e-01 
2334      <output ommitted> 
2335      193   1544 3.15195320957214e-15 1.2e+03 3.70e-08 3.45e-11 
2336      >>> es.logger.disp_header() 
2337      Iterat Nfevals  function value    axis ratio maxstd   minstd 
2338      >>> # es.logger.plot() # will make a plot 
2339   
2340      Example implementing restarts with increasing popsize (IPOP), output 
2341      is not displayed: 
2342   
2343      >>> import cma, numpy as np 
2344      >>> 
2345      >>> # restart with increasing population size (IPOP) 
2346      >>> bestever = cma.BestSolution() 
2347      >>> for lam in 10 * 2**np.arange(8):  # 10, 20, 40, 80, ..., 10 * 2**7 
2348      ...     es = cma.CMAEvolutionStrategy('6 - 8 * np.random.rand(9)',  # 9-D 
2349      ...                                   5,  # initial std sigma0 
2350      ...                                   {'popsize': lam,  # options 
2351      ...                                    'verb_append': bestever.evalsall}) 
2352      ...     logger = cma.CMADataLogger().register(es, append=bestever.evalsall) 
2353      ...     while not es.stop(): 
2354      ...         X = es.ask()    # get list of new solutions 
2355      ...         fit = [cma.fcts.rastrigin(x) for x in X]  # evaluate each solution 
2356      ...         es.tell(X, fit) # besides for termination only the ranking in fit is used 
2357      ... 
2358      ...         # display some output 
2359      ...         logger.add()  # add a "data point" to the log, writing in files 
2360      ...         es.disp()  # uses option verb_disp with default 100 
2361      ... 
2362      ...     print('termination:', es.stop()) 
2363      ...     cma.pprint(es.best.__dict__) 
2364      ... 
2365      ...     bestever.update(es.best) 
2366      ... 
2367      ...     # show a plot 
2368      ...     # logger.plot(); 
2369      ...     if bestever.f < 1e-8:  # global optimum was hit 
2370      ...         break 
2371      <output omitted> 
2372      >>> assert es.result()[1] < 1e-8 
2373   
2374      On the Rastrigin function, usually after five restarts the global optimum 
2375      is located. 
2376   
2377      Using the ``multiprocessing`` module, we can evaluate the function in parallel with a simple 
2378      modification of the example (however multiprocessing seems not always reliable) :: 
2379   
2380          try: 
2381              import multiprocessing as mp 
2382              import cma 
2383              es = cma.CMAEvolutionStrategy(22 * [0.0], 1.0, {'maxiter':10}) 
2384              pool = mp.Pool(es.popsize) 
2385              while not es.stop(): 
2386                  X = es.ask() 
2387                  es.tell(X, pool.map_async(cma.felli, X).get()) # use chunksize parameter as popsize/len(pool)? 
2388                  es.logger.add() 
2389          except ImportError: 
2390              pass 
2391   
2392      The final example shows how to resume: 
2393   
2394      >>> import cma, pickle 
2395      >>> 
2396      >>> es = cma.CMAEvolutionStrategy(12 * [0.1],  # a new instance, 12-D 
2397      ...                               0.5)         # initial std sigma0 
2398      >>> es.optimize(cma.fcts.rosen, iterations=100) 
2399      >>> pickle.dump(es, open('saved-cma-object.pkl', 'wb')) 
2400      >>> print('saved') 
2401      >>> del es  # let's start fresh 
2402      >>> 
2403      >>> es = pickle.load(open('saved-cma-object.pkl', 'rb')) 
2404      >>> print('resumed') 
2405      >>> es.optimize(cma.fcts.rosen, verb_disp=200) 
2406      >>> assert es.result()[2] < 15000 
2407      >>> cma.pprint(es.result()) 
2408   
2409      Details 
2410      ======= 
2411      The following two enhancements are implemented, the latter is turned 
2412      on by default only for very small population size. 
2413   
2414      *Active CMA* is implemented with option ``CMA_active`` and 
2415      conducts an update of the covariance matrix with negative weights. 
2416      The negative update is implemented, such that positive definiteness 
2417      is guarantied. The update is applied after the default update and 
2418      only before the covariance matrix is decomposed, which limits the 
2419      additional computational burden to be at most a factor of three 
2420      (typically smaller). A typical speed up factor (number of 
2421      f-evaluations) is between 1.1 and two. 
2422   
2423      References: Jastrebski and Arnold, CEC 2006, Glasmachers et al, GECCO 2010. 
2424   
2425      *Selective mirroring* is implemented with option ``CMA_mirrors`` 
2426      in the method ``get_mirror()``. Only the method `ask_and_eval()` 
2427      (used by `fmin`) will then sample selectively mirrored vectors. In 
2428      selective mirroring, only the worst solutions are mirrored. With 
2429      the default small number of mirrors, *pairwise selection* (where at 
2430      most one of the two mirrors contribute to the update of the 
2431      distribution mean) is implicitly guarantied under selective 
2432      mirroring and therefore not explicitly implemented. 
2433   
2434      References: Brockhoff et al, PPSN 2010, Auger et al, GECCO 2011. 
2435   
2436      :See: `fmin()`, `CMAOptions`, `plot()`, `ask()`, `tell()`, `ask_and_eval()` 
2437   
2438      """ 
2439      # ____________________________________________________________ 
2440      @property  # read only attribute decorator for a method 
2441 -    def popsize(self): 
2442          """number of samples by default returned by` ask()` 
2443          """ 
2444          return self.sp.popsize 
2445   
2446      # this is not compatible with python2.5: 
2447      #     @popsize.setter 
2448      #     def popsize(self, p): 
2449      #         """popsize cannot be set (this might change in future) 
2450      #         """ 
2451      #         raise _Error("popsize cannot be changed (this might change in future)") 
2452   
2453      # ____________________________________________________________ 
2454      # ____________________________________________________________ 
2455 -    def stop(self, check=True): 
2456          """return a dictionary with the termination status. 
2457          With ``check==False``, the termination conditions are not checked and 
2458          the status might not reflect the current situation. 
2459   
2460          """ 
2461          if (check and self.countiter > 0 and self.opts['termination_callback'] and 
2462                  self.opts['termination_callback'] != str(self.opts['termination_callback'])): 
2463              self.callbackstop = self.opts['termination_callback'](self) 
2464   
2465          return self.stopdict(self if check else None)  # update the stopdict and return a Dict 
2466   
2467      # ____________________________________________________________ 
2468      # ____________________________________________________________ 
2469 -    def __init__(self, x0, sigma0, inopts={}): 
2470          """see class `CMAEvolutionStrategy` 
2471   
2472          """ 
2473          self.inputargs = dict(locals())  # for the record 
2474          del self.inputargs['self']  # otherwise the instance self has a cyclic reference 
2475          self.inopts = inopts 
2476          opts = CMAOptions(inopts).complement()  # CMAOptions() == fmin([],[]) == defaultOptions() 
2477   
2478          if 'noise_handling' in opts and opts.eval('noise_handling'): 
2479              raise ValueError('noise_handling not available with class CMAEvolutionStrategy, use function fmin') 
2480          if 'restarts' in opts and opts.eval('restarts'): 
2481              raise ValueError('restarts not available with class CMAEvolutionStrategy, use function fmin') 
2482   
2483          self.set_x0(x0)  # manage weird shapes, set self.x0 
2484          self.N_pheno = len(self.x0) 
2485   
2486          self.sigma0 = sigma0 
2487          if isinstance(sigma0, str):  # TODO: no real need here (do rather in fmin) 
2488              self.sigma0 = eval(sigma0)  # like '1./N' or 'np.random.rand(1)[0]+1e-2' 
2489          if np.size(self.sigma0) != 1 or np.shape(self.sigma0): 
2490              raise _Error('input argument sigma0 must be (or evaluate to) a scalar') 
2491          self.sigma = self.sigma0  # goes to inialize 
2492   
2493          # extract/expand options 
2494          N = self.N_pheno 
2495          assert isinstance(opts['fixed_variables'], (basestring, dict)) \ 
2496              or opts['fixed_variables'] is None 
2497          # TODO: in case of a string we need to eval the fixed_variables 
2498          if isinstance(opts['fixed_variables'], dict): 
2499              N = self.N_pheno - len(opts['fixed_variables']) 
2500          opts.evalall(locals())  # using only N 
2501          self.opts = opts 
2502   
2503          self.randn = opts['randn'] 
2504          self.gp = GenoPheno(self.N_pheno, opts['scaling_of_variables'], opts['typical_x'], 
2505              opts['fixed_variables'], opts['transformation']) 
2506          self.boundary_handler = opts.eval('boundary_handling')(opts.eval('bounds')) 
2507          if not self.boundary_handler.has_bounds(): 
2508              self.boundary_handler = BoundNone()  # just a little faster and well defined 
2509          elif not self.boundary_handler.is_in_bounds(self.gp.pheno(self.x0)): 
2510              if opts['verbose'] >= 0: 
2511                  print('WARNING: initial solution is out of the domain boundaries:') 
2512                  print('  x0   = ' + str(self.gp.pheno(self.x0))) 
2513                  print('  ldom = ' + str(self.boundary_handler.bounds[0])) 
2514                  print('  udom = ' + str(self.boundary_handler.bounds[1])) 
2515   
2516          # self.mean = array(self.x0, copy=True) 
2517          tmp,  self.gp.tf_geno = self.gp.tf_geno, lambda x: x  # a hack to avoid an exception if tf_geno is None 
2518          self.mean = self.gp.geno(self.x0, copy_always=True) 
2519          self.gp.tf_geno = tmp 
2520          # without copy_always interface: 
2521          # self.mean = self.gp.geno(array(self.x0, copy=True), copy_if_changed=False) 
2522          self.N = len(self.mean) 
2523          assert N == self.N 
2524          self.fmean = np.NaN  # TODO name should change? prints nan in output files (OK with matlab&octave) 
2525          self.fmean_noise_free = 0.  # for output only 
2526   
2527          self.adapt_sigma = opts['AdaptSigma'] 
2528          if self.adapt_sigma is False: 
2529              self.adapt_sigma = CMAAdaptSigmaNone 
2530          self.adapt_sigma = self.adapt_sigma()  # class instance 
2531   
2532          self.sp = _CMAParameters(N, opts) 
2533          self.sp0 = self.sp  # looks useless, as it is not a copy 
2534   
2535          # initialization of state variables 
2536          self.countiter = 0 
2537          self.countevals = max((0, opts['verb_append'])) \ 
2538              if not isinstance(opts['verb_append'], bool) else 0 
2539          self.pc = np.zeros(N) 
2540   
2541          self.sigma_vec = np.ones(N) if np.isfinite(self.sp.dampsvec) else 1 
2542          stds = np.ones(N) 
2543          if self.opts['CMA_teststds'] is not None and np.all(self.opts['CMA_teststds']):  # also 0 would not make sense 
2544              stds = array(self.opts['CMA_teststds']) 
2545              if np.size(stds) != N: 
2546                  raise _Error('CMA_teststds option must have dimension = ' + str(N)) 
2547          if self.opts['CMA_diagonal']:  # is True or > 0 
2548              # linear time and space complexity 
2549              self.B = array(1)  # works fine with np.dot(self.B, anything) and self.B.T 
2550              self.C = stds**2  # TODO: remove this!? 
2551              self.dC = self.C 
2552          else: 
2553              self.B = np.eye(N)  # identity(N), do not from matlib import *, as eye is a matrix there 
2554              # prevent equal eigenvals, a hack for np.linalg: 
2555              self.C = np.diag(stds**2 * exp(1e-6 * (np.random.rand(N) - 0.5))) 
2556              self.dC = np.diag(self.C).copy() 
2557              self.Yneg = np.zeros((N, N)) 
2558          self.D = stds 
2559   
2560          # self.gp.pheno adds fixed variables 
2561          relative_stds = ((self.gp.pheno(self.mean + self.sigma * self.sigma_vec * self.D) 
2562                            - self.gp.pheno(self.mean - self.sigma * self.sigma_vec * self.D)) / 2.0 
2563                           / (self.boundary_handler.get_bounds('upper', self.N_pheno) 
2564                              - self.boundary_handler.get_bounds('lower', self.N_pheno))) 
2565          if np.any(relative_stds > 1): 
2566              raise ValueError('initial standard deviations larger than the bounded domain size in variables ' 
2567                           + str(np.where(relative_stds > 1)[0])) 
2568          self.flgtelldone = True 
2569          self.itereigenupdated = self.countiter 
2570          self.noiseS = 0  # noise "signal" 
2571          self.hsiglist = [] 
2572   
2573          if not opts['seed']: 
2574              np.random.seed() 
2575              six_decimals = (time.time() - 1e6 * (time.time() // 1e6)) 
2576              opts['seed'] = 1e5 * np.random.rand() + six_decimals + 1e5 * (time.time() % 1) 
2577          opts['seed'] = int(opts['seed']) 
2578          np.random.seed(opts['seed']) 
2579   
2580          self.sent_solutions = CMASolutionDict() 
2581          self.archive = CMASolutionDict() 
2582          self.best = BestSolution() 
2583   
2584          out = {}  # TODO: obsolete, replaced by method results()? 
2585          out['best'] = self.best 
2586          # out['hsigcount'] = 0 
2587          out['termination'] = {} 
2588          self.out = out 
2589   
2590          self.const = _BlancClass() 
2591          self.const.chiN = N**0.5 * (1 - 1. / (4.*N) + 1. / (21.*N**2))  # expectation of norm(randn(N,1)) 
2592   
2593          self.logger = CMADataLogger(opts['verb_filenameprefix'], modulo=opts['verb_log']).register(self) 
2594   
2595          # attribute for stopping criteria in function stop 
2596          self.stopdict = CMAStopDict() 
2597          self.callbackstop = 0 
2598   
2599          self.fit = _BlancClass() 
2600          self.fit.fit = []  # not really necessary 
2601          self.fit.hist = []  # short history of best 
2602          self.fit.histbest = []  # long history of best 
2603          self.fit.histmedian = []  # long history of median 
2604   
2605          self.more_to_write = []  # [1, 1, 1, 1]  #  N*[1]  # needed when writing takes place before setting 
2606   
2607          # say hello 
2608          if opts['verb_disp'] > 0 and opts['verbose'] >= 0: 
2609              sweighted = '_w' if self.sp.mu > 1 else '' 
2610              smirr = 'mirr%d' % (self.sp.lam_mirr) if self.sp.lam_mirr else '' 
2611              print('(%d' % (self.sp.mu) + sweighted + ',%d' % (self.sp.popsize) + smirr + 
2612                    ')-' + ('a' if opts['CMA_active'] else '') + 'CMA-ES' + 
2613                    ' (mu_w=%2.1f,w_1=%d%%)' % (self.sp.mueff, int(100 * self.sp.weights[0])) + 
2614                    ' in dimension %d (seed=%d, %s)' % (N, opts['seed'], time.asctime()))  # + func.__name__ 
2615              if opts['CMA_diagonal'] and self.sp.CMA_on: 
2616                  s = '' 
2617                  if opts['CMA_diagonal'] is not True: 
2618                      s = ' for ' 
2619                      if opts['CMA_diagonal'] < np.inf: 
2620                          s += str(int(opts['CMA_diagonal'])) 
2621                      else: 
2622                          s += str(np.floor(opts['CMA_diagonal'])) 
2623                      s += ' iterations' 
2624                      s += ' (1/ccov=' + str(round(1. / (self.sp.c1 + self.sp.cmu))) + ')' 
2625                  print('   Covariance matrix is diagonal' + s) 
2626   
2627 -    def set_x0(self, x0): 
2628          if x0 == str(x0): 
2629              x0 = eval(x0) 
2630          self.x0 = array(x0)  # should not have column or row, is just 1-D 
2631          if self.x0.ndim == 2: 
2632              if self.opts.eval('verbose') >= 0: 
2633                  print('WARNING: input x0 should be a list or 1-D array, trying to flatten ' + 
2634                          str(self.x0.shape) + '-array') 
2635              if self.x0.shape[0] == 1: 
2636                  self.x0 = self.x0[0] 
2637              elif self.x0.shape[1] == 1: 
2638                  self.x0 = array([x[0] for x in self.x0]) 
2639          if self.x0.ndim != 1: 
2640              raise _Error('x0 must be 1-D array') 
2641          if len(self.x0) <= 1: 
2642              raise _Error('optimization in 1-D is not supported (code was never tested)') 
2643          self.x0.resize(self.x0.shape[0])  # 1-D array, not really necessary?! 
2644   
2645      # ____________________________________________________________ 
2646      # ____________________________________________________________ 
2647 -    def ask(self, number=None, xmean=None, sigma_fac=1): 
2648          """get new candidate solutions, sampled from a multi-variate 
2649          normal distribution and transformed to f-representation 
2650          (phenotype) to be evaluated. 
2651   
2652          Arguments 
2653          --------- 
2654              `number` 
2655                  number of returned solutions, by default the 
2656                  population size ``popsize`` (AKA ``lambda``). 
2657              `xmean` 
2658                  distribution mean 
2659              `sigma_fac` 
2660                  multiplier for internal sample width (standard 
2661                  deviation) 
2662   
2663          Return 
2664          ------ 
2665          A list of N-dimensional candidate solutions to be evaluated 
2666   
2667          Example 
2668          ------- 
2669          >>> import cma 
2670          >>> es = cma.CMAEvolutionStrategy([0,0,0,0], 0.3) 
2671          >>> while not es.stop() and es.best.f > 1e-6:  # my_desired_target_f_value 
2672          ...     X = es.ask()  # get list of new solutions 
2673          ...     fit = [cma.fcts.rosen(x) for x in X]  # call function rosen with each solution 
2674          ...     es.tell(X, fit)  # feed values 
2675   
2676          :See: `ask_and_eval`, `ask_geno`, `tell` 
2677   
2678          """ 
2679          pop_geno = self.ask_geno(number, xmean, sigma_fac) 
2680   
2681   
2682          # N,lambda=20,200: overall CPU 7s vs 5s == 40% overhead, even without bounds! 
2683          #                  new data: 11.5s vs 9.5s == 20% 
2684          # TODO: check here, whether this is necessary? 
2685          # return [self.gp.pheno(x, copy=False, into_bounds=self.boundary_handler.repair) for x in pop]  # probably fine 
2686          # return [Solution(self.gp.pheno(x, copy=False), copy=False) for x in pop]  # here comes the memory leak, now solved 
2687          # pop_pheno = [Solution(self.gp.pheno(x, copy=False), copy=False).repair(self.gp.bounds) for x in pop_geno] 
2688          pop_pheno = [self.gp.pheno(x, copy=True, into_bounds=self.boundary_handler.repair) for x in pop_geno] 
2689   
2690          # insert solutions, this could also (better?) be done in self.gp.pheno 
2691          for i in rglen((pop_geno)): 
2692              self.sent_solutions.insert(pop_pheno[i], geno=pop_geno[i], iteration=self.countiter) 
2693          return pop_pheno 
2694   
2695      # ____________________________________________________________ 
2696      # ____________________________________________________________ 
2697 -    def ask_geno(self, number=None, xmean=None, sigma_fac=1): 
2698          """get new candidate solutions in genotyp, sampled from a 
2699          multi-variate normal distribution. 
2700   
2701          Arguments are 
2702              `number` 
2703                  number of returned solutions, by default the 
2704                  population size `popsize` (AKA lambda). 
2705              `xmean` 
2706                  distribution mean 
2707              `sigma_fac` 
2708                  multiplier for internal sample width (standard 
2709                  deviation) 
2710   
2711          `ask_geno` returns a list of N-dimensional candidate solutions 
2712          in genotyp representation and is called by `ask`. 
2713   
2714          :See: `ask`, `ask_and_eval` 
2715   
2716          """ 
2717   
2718          if number is None or number < 1: 
2719              number = self.sp.popsize 
2720          if xmean is None: 
2721              xmean = self.mean 
2722          else: 
2723              try: 
2724                  xmean = self.archive[xmean]['geno'] 
2725                  # noise handling after call of tell 
2726              except KeyError: 
2727                  try: 
2728                      xmean = self.sent_solutions[xmean]['geno'] 
2729                      # noise handling before calling tell 
2730                  except KeyError: 
2731                      pass 
2732   
2733          if self.countiter == 0: 
2734              self.tic = time.clock()  # backward compatible 
2735              self.elapsed_time = ElapsedTime() 
2736   
2737          sigma = sigma_fac * self.sigma 
2738   
2739          # update parameters for sampling the distribution 
2740          #        fac  0      1      10 
2741          # 150-D cigar: 
2742          #           50749  50464   50787 
2743          # 200-D elli:               == 6.9 
2744          #                  99900   101160 
2745          #                 100995   103275 == 2% loss 
2746          # 100-D elli:               == 6.9 
2747          #                 363052   369325  < 2% loss 
2748          #                 365075   365755 
2749   
2750          # update distribution 
2751          if self.sp.CMA_on and ( 
2752                  (self.opts['updatecovwait'] is None and 
2753                   self.countiter >= 
2754                       self.itereigenupdated + 1. / (self.sp.c1 + self.sp.cmu) / self.N / 10 
2755                   ) or 
2756                  (self.opts['updatecovwait'] is not None and 
2757                   self.countiter > self.itereigenupdated + self.opts['updatecovwait'] 
2758                   ) or 
2759                  (self.sp.neg.cmuexp * (self.countiter - self.itereigenupdated) > 0.5 
2760                  )  # TODO (minor): not sure whether this is "the right" criterion 
2761              ): 
2762              self.updateBD() 
2763          # sample distribution 
2764          if self.flgtelldone:  # could be done in tell()!? 
2765              self.flgtelldone = False 
2766              self.ary = [] 
2767   
2768          # each row is a solution 
2769          arz = self.randn((number, self.N)) 
2770          # zzzzzzzzzzzzzzzzzzzzzzzzzzz 
2771          if self.opts['CMA_sample_on_sphere_surface']:  # normalize the length to chiN 
2772              for i in rglen((arz)): 
2773                  ss = sum(arz[i]**2) 
2774                  if 1 < 3 or ss > self.N + 10.1: 
2775                      arz[i] *= (self.N**0.5 if self.opts['CSA_squared'] else self.const.chiN) / ss**0.5 
2776              # or to average 
2777              # arz *= 1 * self.const.chiN / np.mean([sum(z**2)**0.5 for z in arz]) 
2778   
2779          # fac = np.mean(sum(arz**2, 1)**0.5) 
2780          # print fac 
2781          # arz *= self.const.chiN / fac 
2782          ary = self.sigma_vec * np.dot(self.B, (self.D * arz).T).T 
2783          if number > 2 and self.countiter > 2: 
2784              if (isinstance(self.adapt_sigma, CMAAdaptSigmaTPA) or 
2785                  self.opts['mean_shift_line_samples'] or 
2786                  self.opts['pc_line_samples']): 
2787                  ys = [] 
2788                  if self.opts['pc_line_samples']: 
2789                      ys.append(self.pc[:])  # now TPA is with pc_line_samples 
2790                  if self.opts['mean_shift_line_samples']: 
2791                      ys.append(self.mean - self.mean_old) 
2792                  if not len(ys): 
2793                      ys.append(self.mean - self.mean_old) 
2794                  # assign a mirrored pair from each element of ys into ary 
2795                  for i, y in enumerate(ys): 
2796                      if len(arz) > 2 * i + 1:  # at least two more samples 
2797                          assert y is not self.pc 
2798                          y *= sum(self.randn(self.N)**2)**0.5 / self.mahalanobisNorm(y) 
2799                          # TODO: rescale y depending on some parameter? 
2800                          ary[2*i] = y / self.sigma 
2801                          ary[2*i + 1] = y / -self.sigma 
2802                      else: 
2803                          _print_warning('line samples omitted due to small popsize', 
2804                              method_name='ask_geno', iteration=self.countiter) 
2805   
2806          pop = xmean + sigma * ary 
2807          self.evaluations_per_f_value = 1 
2808          self.ary = ary  # ask_geno is called recursively in CMAAdaptSigmaTPA 
2809          if number == self.sp.popsize: 
2810              self.arz = arz  # is never used 
2811          return pop 
2812   
2813 -    def get_mirror(self, x, preserve_length=False): 
2814          """return ``pheno(self.mean - (geno(x) - self.mean))``. 
2815   
2816          >>> import cma 
2817          >>> es = cma.CMAEvolutionStrategy(cma.np.random.randn(3), 1) 
2818          >>> x = cma.np.random.randn(3) 
2819          >>> assert cma.Mh.vequals_approximately(es.mean - (x - es.mean), es.get_mirror(x, preserve_length=True)) 
2820          >>> x = es.ask(1)[0] 
2821          >>> vals = (es.get_mirror(x) - es.mean) / (x - es.mean) 
2822          >>> assert cma.Mh.equals_approximately(sum(vals), len(vals) * vals[0]) 
2823   
2824          TODO: this implementation is yet experimental. 
2825   
2826          Selectively mirrored sampling improves to a moderate extend but 
2827          overadditively with active CMA for quite understandable reasons. 
2828   
2829          Optimal number of mirrors are suprisingly small: 1,2,3 for maxlam=7,13,20 
2830          however note that 3,6,10 are the respective maximal possible mirrors that 
2831          must be clearly suboptimal. 
2832   
2833          """ 
2834          try: 
2835              dx = self.sent_solutions[x]['geno'] - self.mean 
2836          except:  # can only happen with injected solutions?! 
2837              dx = self.gp.geno(x, from_bounds=self.boundary_handler.inverse, copy_if_changed=True) - self.mean 
2838   
2839          if not preserve_length: 
2840              dx *= sum(self.randn(self.N)**2)**0.5 / self.mahalanobisNorm(dx) 
2841          x = self.mean - dx 
2842          y = self.gp.pheno(x, into_bounds=self.boundary_handler.repair) 
2843          # old measure: costs 25% in CPU performance with N,lambda=20,200 
2844          self.sent_solutions.insert(y, geno=x, iteration=self.countiter) 
2845          return y 
2846   
2847 -    def mirror_penalized(self, f_values, idx): 
2848          """obsolete and subject to removal (TODO), 
2849          return modified f-values such that for each mirror one becomes worst. 
2850   
2851          This function is useless when selective mirroring is applied with no 
2852          more than (lambda-mu)/2 solutions. 
2853   
2854          Mirrors are leading and trailing values in ``f_values``. 
2855   
2856          """ 
2857          assert len(f_values) >= 2 * len(idx) 
2858          m = np.max(np.abs(f_values)) 
2859          for i in len(idx): 
2860              if f_values[idx[i]] > f_values[-1 - i]: 
2861                  f_values[idx[i]] += m 
2862              else: 
2863                  f_values[-1 - i] += m 
2864          return f_values 
2865   
2866 -    def mirror_idx_cov(self, f_values, idx1):  # will most likely be removed 
2867          """obsolete and subject to removal (TODO), 
2868          return indices for negative ("active") update of the covariance matrix 
2869          assuming that ``f_values[idx1[i]]`` and ``f_values[-1-i]`` are 
2870          the corresponding mirrored values 
2871   
2872          computes the index of the worse solution sorted by the f-value of the 
2873          better solution. 
2874   
2875          TODO: when the actual mirror was rejected, it is better 
2876          to return idx1 instead of idx2. 
2877   
2878          Remark: this function might not be necessary at all: if the worst solution 
2879          is the best mirrored, the covariance matrix updates cancel (cave: weights 
2880          and learning rates), which seems what is desirable. If the mirror is bad, 
2881          as strong negative update is made, again what is desirable. 
2882          And the fitness--step-length correlation is in part addressed by 
2883          using flat weights. 
2884   
2885          """ 
2886          idx2 = np.arange(len(f_values) - 1, len(f_values) - 1 - len(idx1), -1) 
2887          f = [] 
2888          for i in rglen((idx1)): 
2889              f.append(min((f_values[idx1[i]], f_values[idx2[i]]))) 
2890              # idx.append(idx1[i] if f_values[idx1[i]] > f_values[idx2[i]] else idx2[i]) 
2891          return idx2[np.argsort(f)][-1::-1] 
2892   
2893 -    def eval_mean(self, func, args=()): 
2894          """evaluate the distribution mean, this is not (yet) effective 
2895          in terms of termination or display""" 
2896          self.fmean = func(self.mean, *args) 
2897   
2898      # ____________________________________________________________ 
2899      # ____________________________________________________________ 
2900      # 
2901 -    def ask_and_eval(self, func, args=(), number=None, xmean=None, sigma_fac=1, 
2902                       evaluations=1, aggregation=np.median, kappa=1): 
2903          """samples `number` solutions and evaluates them on `func`, where 
2904          each solution `s` is resampled until ``self.is_feasible(s, func(s)) is True``. 
2905   
2906          Arguments 
2907          --------- 
2908              `func` 
2909                  objective function 
2910              `args` 
2911                  additional parameters for `func` 
2912              `number` 
2913                  number of solutions to be sampled, by default 
2914                  population size ``popsize`` (AKA lambda) 
2915              `xmean` 
2916                  mean for sampling the solutions, by default ``self.mean``. 
2917              `sigma_fac` 
2918                  multiplier for sampling width, standard deviation, for example 
2919                  to get a small perturbation of solution `xmean` 
2920              `evaluations` 
2921                  number of evaluations for each sampled solution 
2922              `aggregation` 
2923                  function that aggregates `evaluations` values to 
2924                  as single value. 
2925              `kappa` 
2926                  multiplier used for the evaluation of the solutions, in 
2927                  that ``func(m + kappa*(x - m))`` is the f-value for x. 
2928   
2929          Return 
2930          ------ 
2931          ``(X, fit)``, where 
2932              X -- list of solutions 
2933              fit -- list of respective function values 
2934   
2935          Details 
2936          ------- 
2937          While ``not self.is_feasible(x, func(x))``new solutions are sampled. By 
2938          default ``self.is_feasible == cma.feasible == lambda x, f: f not in (None, np.NaN)``. 
2939          The argument to `func` can be freely modified within `func`. 
2940   
2941          Depending on the ``CMA_mirrors`` option, some solutions are not sampled 
2942          independently but as mirrors of other bad solutions. This is a simple 
2943          derandomization that can save 10-30% of the evaluations in particular 
2944          with small populations, for example on the cigar function. 
2945   
2946          Example 
2947          ------- 
2948          >>> import cma 
2949          >>> x0, sigma0 = 8*[10], 1  # 8-D 
2950          >>> es = cma.CMAEvolutionStrategy(x0, sigma0) 
2951          >>> while not es.stop(): 
2952          ...     X, fit = es.ask_and_eval(cma.fcts.elli)  # handles NaN with resampling 
2953          ...     es.tell(X, fit)  # pass on fitness values 
2954          ...     es.disp(20) # print every 20-th iteration 
2955          >>> print('terminated on ' + str(es.stop())) 
2956          <output omitted> 
2957   
2958          A single iteration step can be expressed in one line, such that 
2959          an entire optimization after initialization becomes 
2960          :: 
2961   
2962              while not es.stop(): 
2963                  es.tell(*es.ask_and_eval(cma.fcts.elli)) 
2964   
2965          """ 
2966          # initialize 
2967          popsize = self.sp.popsize 
2968          if number is not None: 
2969              popsize = number 
2970          selective_mirroring = True 
2971          nmirrors = self.sp.lam_mirr 
2972          if popsize != self.sp.popsize: 
2973              nmirrors = Mh.sround(popsize * self.sp.lam_mirr / self.sp.popsize) 
2974              # TODO: now selective mirroring might be impaired 
2975          assert nmirrors <= popsize // 2 
2976          self.mirrors_idx = np.arange(nmirrors)  # might never be used 
2977          self.mirrors_rejected_idx = []  # might never be used 
2978          if xmean is None: 
2979              xmean = self.mean 
2980          is_feasible = self.opts['is_feasible'] 
2981   
2982          # do the work 
2983          fit = []  # or np.NaN * np.empty(number) 
2984          X_first = self.ask(popsize) 
2985          X = [] 
2986          for k in xrange(int(popsize)): 
2987              x, f = X_first.pop(0), None 
2988              nreject = -1 
2989              while nreject < 0 or not is_feasible(x, f):  # rejection sampling 
2990                  nreject += 1 
2991                  if nreject:  # resample 
2992                      x = self.ask(1, xmean, sigma_fac)[0] 
2993                  elif k >= popsize - nmirrors:  # mirrored sample 
2994                      if k == popsize - nmirrors and selective_mirroring: 
2995                          self.mirrors_idx = np.argsort(fit)[-1:-1 - nmirrors:-1] 
2996                      x = self.get_mirror(X[self.mirrors_idx[popsize - 1 - k]]) 
2997                  if nreject == 1 and k >= popsize - nmirrors: 
2998                      self.mirrors_rejected_idx.append(k) 
2999   
3000                  # contraints handling test hardwired ccccccccccc 
3001                  length_normalizer = 1 
3002                  # zzzzzzzzzzzzzzzzzzzzzzzzz 
3003                  f = func(x, *args) if kappa == 1 else func(xmean + kappa * length_normalizer * (x - xmean), *args) 
3004                  if is_feasible(x, f) and evaluations > 1: 
3005                      f = aggregation([f] + [(func(x, *args) if kappa == 1 else func(xmean + kappa * length_normalizer * (x - xmean), *args)) for _i in xrange(int(evaluations - 1))]) 
3006                  if nreject + 1 % 1000 == 0: 
3007                      print('  %d solutions rejected (f-value NaN or None) at iteration %d' % 
3008                            (nreject, self.countiter)) 
3009              fit.append(f) 
3010              X.append(x) 
3011          self.evaluations_per_f_value = int(evaluations) 
3012          return X, fit 
3013   
3014   
3015      # ____________________________________________________________ 
3016 -    def tell(self, solutions, function_values, check_points=None, copy=False): 
3017          """pass objective function values to prepare for next 
3018          iteration. This core procedure of the CMA-ES algorithm updates 
3019          all state variables, in particular the two evolution paths, the 
3020          distribution mean, the covariance matrix and a step-size. 
3021   
3022          Arguments 
3023          --------- 
3024              `solutions` 
3025                  list or array of candidate solution points (of 
3026                  type `numpy.ndarray`), most presumably before 
3027                  delivered by method `ask()` or `ask_and_eval()`. 
3028              `function_values` 
3029                  list or array of objective function values 
3030                  corresponding to the respective points. Beside for termination 
3031                  decisions, only the ranking of values in `function_values` 
3032                  is used. 
3033              `check_points` 
3034                  If ``check_points is None``, only solutions that are not generated 
3035                  by `ask()` are possibly clipped (recommended). ``False`` does not clip 
3036                  any solution (not recommended). 
3037                  If ``True``, clips solutions that realize long steps (i.e. also 
3038                  those that are unlikely to be generated with `ask()`). `check_points` 
3039                  can be a list of indices to be checked in solutions. 
3040              `copy` 
3041                  ``solutions`` can be modified in this routine, if ``copy is False`` 
3042              `AdaptSigma`: 
3043                  sigma adaptation class like ``CMAAdaptSigmaCSA``, with an adhoc interface 
3044                  very specific to the ``CMAEvolutionStrategy.tell`` method 
3045                  (this interface might change in future). Overwrites `self.AdaptSigma`. 
3046   
3047          Details 
3048          ------- 
3049          `tell()` updates the parameters of the multivariate 
3050          normal search distribution, namely covariance matrix and 
3051          step-size and updates also the attributes `countiter` and 
3052          `countevals`. To check the points for consistency is quadratic 
3053          in the dimension (like sampling points). 
3054   
3055          Bugs 
3056          ---- 
3057          The effect of changing the solutions delivered by `ask()` depends on whether 
3058          boundary handling is applied. With boundary handling, modifications are 
3059          disregarded. This is necessary to apply the default boundary handling that 
3060          uses unrepaired solutions but might change in future. 
3061   
3062          Example 
3063          ------- 
3064          :: 
3065   
3066              import cma 
3067              func = cma.fcts.elli  # choose objective function 
3068              es = cma.CMAEvolutionStrategy(cma.np.random.rand(10), 1) 
3069              while not es.stop(): 
3070                 X = es.ask() 
3071                 es.tell(X, [func(x) for x in X]) 
3072              es.result()  # where the result can be found 
3073   
3074          :See: class `CMAEvolutionStrategy`, `ask()`, `ask_and_eval()`, `fmin()` 
3075   
3076          """ 
3077          if self.flgtelldone: 
3078              raise _Error('tell should only be called once per iteration') 
3079   
3080          lam = len(solutions) 
3081          if lam != array(function_values).shape[0]: 
3082              raise _Error('for each candidate solution ' 
3083                          + 'a function value must be provided') 
3084          if lam + self.sp.lam_mirr < 3: 
3085              raise _Error('population size ' + str(lam) + ' is too small when option CMA_mirrors * popsize < 0.5') 
3086   
3087          if not np.isscalar(function_values[0]): 
3088              if np.isscalar(function_values[0][0]): 
3089                  if self.countiter <= 1: 
3090                      print('WARNING: function values are not a list of scalars (further warnings are suppressed)') 
3091                  function_values = [val[0] for val in function_values] 
3092              else: 
3093                  raise _Error('objective function values must be a list of scalars') 
3094   
3095   
3096          # ## prepare 
3097          N = self.N 
3098          sp = self.sp 
3099          if lam < sp.mu:  # rather decrease cmean instead of having mu > lambda//2 
3100              raise _Error('not enough solutions passed to function tell (mu>lambda)') 
3101   
3102          self.countiter += 1  # >= 1 now 
3103          self.countevals += sp.popsize * self.evaluations_per_f_value 
3104          self.best.update(solutions, self.sent_solutions, function_values, self.countevals) 
3105   
3106          flgseparable = self.opts['CMA_diagonal'] is True \ 
3107                         or self.countiter <= self.opts['CMA_diagonal'] 
3108          if not flgseparable and len(self.C.shape) == 1:  # C was diagonal ie 1-D 
3109              # enter non-separable phase (no easy return from here) 
3110              self.B = np.eye(N)  # identity(N) 
3111              self.C = np.diag(self.C) 
3112              idx = np.argsort(self.D) 
3113              self.D = self.D[idx] 
3114              self.B = self.B[:, idx] 
3115              self.Yneg = np.zeros((N, N)) 
3116   
3117          # ## manage fitness 
3118          fit = self.fit  # make short cut 
3119   
3120          # CPU for N,lam=20,200: this takes 10s vs 7s 
3121          fit.bndpen = self.boundary_handler.update(function_values, self)(solutions, self.sent_solutions, self.gp) 
3122          # for testing: 
3123          # fit.bndpen = self.boundary_handler.update(function_values, self)([s.unrepaired for s in solutions]) 
3124          fit.idx = np.argsort(array(fit.bndpen) + array(function_values)) 
3125          fit.fit = array(function_values, copy=False)[fit.idx] 
3126   
3127          # update output data TODO: this is obsolete!? However: need communicate current best x-value? 
3128          # old: out['recent_x'] = self.gp.pheno(pop[0]) 
3129          self.out['recent_x'] = array(solutions[fit.idx[0]])  # TODO: change in a data structure(?) and use current as identify 
3130          self.out['recent_f'] = fit.fit[0] 
3131   
3132          # fitness histories 
3133          fit.hist.insert(0, fit.fit[0]) 
3134          # if len(self.fit.histbest) < 120+30*N/sp.popsize or  # does not help, as tablet in the beginning is the critical counter-case 
3135          if ((self.countiter % 5) == 0):  # 20 percent of 1e5 gen. 
3136              fit.histbest.insert(0, fit.fit[0]) 
3137              fit.histmedian.insert(0, np.median(fit.fit) if len(fit.fit) < 21 
3138                                      else fit.fit[self.popsize // 2]) 
3139          if len(fit.histbest) > 2e4:  # 10 + 30*N/sp.popsize: 
3140              fit.histbest.pop() 
3141              fit.histmedian.pop() 
3142          if len(fit.hist) > 10 + 30 * N / sp.popsize: 
3143              fit.hist.pop() 
3144   
3145          # TODO: clean up inconsistency when an unrepaired solution is available and used 
3146          # now get the genotypes 
3147          pop = []  # create pop from input argument solutions 
3148          for k, s in enumerate(solutions):  # use phenotype before Solution.repair() 
3149              if 1 < 3: 
3150                  pop += [self.gp.geno(s, 
3151                                       from_bounds=self.boundary_handler.inverse, 
3152                                       repair=(self.repair_genotype if check_points not in (False, 0, [], ()) else None), 
3153                                       archive=self.sent_solutions)]  # takes genotype from sent_solutions, if available 
3154                  try: 
3155                      self.archive.insert(s, value=self.sent_solutions.pop(s), fitness=function_values[k]) 
3156                      # self.sent_solutions.pop(s) 
3157                  except KeyError: 
3158                      pass 
3159          try: 
3160              moldold = self.mean_old 
3161          except: 
3162              pass 
3163          self.mean_old = self.mean 
3164          mold = self.mean_old  # just an alias 
3165   
3166          # check and normalize each x - m 
3167          # check_points is a flag (None is default: check non-known solutions) or an index list 
3168          # should also a number possible (first check_points points)? 
3169          if check_points not in (None, False, 0, [], ()):  # useful in case of injected solutions and/or adaptive encoding, however is automatic with use_sent_solutions 
3170              try: 
3171                  if len(check_points): 
3172                      idx = check_points 
3173              except: 
3174                  idx = xrange(sp.popsize) 
3175   
3176              for k in idx: 
3177                  self.repair_genotype(pop[k]) 
3178   
3179          # only arrays can be multiple indexed 
3180          pop = array(pop, copy=False) 
3181   
3182          # sort pop 
3183          pop = pop[fit.idx] 
3184   
3185          if self.opts['CMA_elitist'] and self.best.f < fit.fit[0]: 
3186              if self.best.x_geno is not None: 
3187                  xp = [self.best.x_geno] 
3188                  # xp = [self.best.xdict['geno']] 
3189                  # xp = [self.gp.geno(self.best.x[:])]  # TODO: remove 
3190                  # print self.mahalanobisNorm(xp[0]-self.mean) 
3191              else: 
3192                  xp = [self.gp.geno(array(self.best.x, copy=True), self.boundary_handler.inverse, copy_if_changed=False)] 
3193                  print('genotype for elitist not found') 
3194              self.clip_or_fit_solutions(xp, [0]) 
3195              pop = array([xp[0]] + list(pop)) 
3196   
3197          # compute new mean 
3198          self.mean = mold + self.sp.cmean * \ 
3199                      (sum(sp.weights * pop[0:sp.mu].T, 1) - mold) 
3200   
3201   
3202          # check Delta m (this is not default, but could become at some point) 
3203          # CAVE: upper_length=sqrt(2)+2 is too restrictive, test upper_length = sqrt(2*N) thoroughly. 
3204          # replaced by repair_geno? 
3205          # simple test case injecting self.mean: 
3206          # self.mean = 1e-4 * self.sigma * np.random.randn(N) 
3207          if 1 < 3: 
3208              cmean = self.sp.cmean 
3209   
3210          # zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz 
3211          # get learning rate constants 
3212          cc, c1, cmu = sp.cc, sp.c1, sp.cmu 
3213          if flgseparable: 
3214              cc, c1, cmu = sp.cc_sep, sp.c1_sep, sp.cmu_sep 
3215   
3216          # now the real work can start 
3217   
3218          hsig = self.adapt_sigma.hsig(self) # ps update must be done here in separable case 
3219   
3220          # hsig = sum(self.ps**2) / self.N < 2 + 4./(N+1) 
3221          # adjust missing variance due to hsig, in 4-D with damps=1e99 and sig0 small 
3222          #       hsig leads to premature convergence of C otherwise 
3223          # hsiga = (1-hsig**2) * c1 * cc * (2-cc)  # to be removed in future 
3224          c1a = c1 - (1 - hsig**2) * c1 * cc * (2 - cc)  # adjust for variance loss 
3225   
3226          self.pc = (1 - cc) * self.pc + \ 
3227                    hsig * (sqrt(cc * (2 - cc) * sp.mueff) / self.sigma / cmean) * \ 
3228                    (self.mean - mold) / self.sigma_vec 
3229   
3230          # covariance matrix adaptation/udpate 
3231          if sp.CMA_on: 
3232              # assert sp.c1 + sp.cmu < sp.mueff / N  # ?? 
3233              assert c1 + cmu <= 1 
3234   
3235              # default full matrix case 
3236              if not flgseparable: 
3237                  Z = (pop[0:sp.mu] - mold) / (self.sigma * self.sigma_vec) 
3238                  Z = dot((cmu * sp.weights) * Z.T, Z)  # learning rate integrated 
3239                  if self.sp.neg.cmuexp: 
3240                      tmp = (pop[-sp.neg.mu:] - mold) / (self.sigma * self.sigma_vec) 
3241                      if 1 < 3:  # normalize to constant length (seems preferable in several aspects) 
3242                          # print(tmp.shape) 
3243                          for i in range(tmp.shape[0]): 
3244                              tmp[i, :] *= N**0.5 / self.mahalanobisNorm(tmp[i, :]) / (self.sigma * self.sigma_vec) 
3245                          # print(tmp.shape) 
3246                      self.Yneg *= 1 - self.sp.neg.cmuexp  # for some reason necessary? 
3247                      self.Yneg += dot(sp.neg.weights * tmp.T, tmp) - self.C 
3248                      # self.update_exponential(dot(sp.neg.weights * tmp.T, tmp) - 1 * self.C, -1*self.sp.neg.cmuexp) 
3249   
3250                  self.C *= 1 - c1a - cmu 
3251                  self.C += np.outer(c1 * self.pc, self.pc) + Z 
3252                  self.dC = np.diag(self.C).copy()  # for output and termination checking 
3253   
3254              else:  # separable/diagonal linear case 
3255                  assert(c1 + cmu <= 1) 
3256                  Z = np.zeros(N) 
3257                  for k in xrange(sp.mu): 
3258                      z = (pop[k] - mold) / (self.sigma * self.sigma_vec)  # TODO see above 
3259                      Z += sp.weights[k] * z * z  # is 1-D 
3260                  self.C = (1 - c1a - cmu) * self.C + c1 * self.pc * self.pc + cmu * Z 
3261                  # TODO: self.C *= exp(cmuneg * (N - dot(sp.neg.weights,  **2) 
3262                  self.dC = self.C 
3263                  self.D = sqrt(self.C)  # C is a 1-D array, this is why adapt_sigma needs to prepare before 
3264                  self.itereigenupdated = self.countiter 
3265   
3266                  # idx = self.mirror_idx_cov()  # take half of mirrored vectors for negative update 
3267   
3268          # step-size adaptation, adapt sigma 
3269          # in case of TPA, function_values[0] and [1] must reflect samples colinear to xmean - xmean_old 
3270          self.adapt_sigma.update(self, function_values=function_values) 
3271   
3272          if self.sigma * min(self.dC)**0.5 < self.opts['minstd']: 
3273              self.sigma = self.opts['minstd'] / min(self.dC)**0.5 
3274          # g = self.countiter 
3275          # N = self.N 
3276          mindx = eval(self.opts['mindx']) if isinstance(self.opts['mindx'], basestring) else self.opts['mindx'] 
3277          if self.sigma * min(self.D) < mindx:  # TODO: sigma_vec is missing here 
3278              self.sigma = mindx / min(self.D) 
3279   
3280          if self.sigma > 1e9 * self.sigma0: 
3281              alpha = self.sigma / max(self.D) 
3282              self.multiplyC(alpha) 
3283              self.sigma /= alpha**0.5 
3284              self.opts['tolupsigma'] /= alpha**0.5  # to be compared with sigma 
3285   
3286          # TODO increase sigma in case of a plateau? 
3287   
3288          # Uncertainty noise measurement is done on an upper level 
3289   
3290          self.flgtelldone = True 
3291      # end tell() 
3292   
3293 -    def result(self): 
3294          """return ``(xbest, f(xbest), evaluations_xbest, evaluations, iterations, pheno(xmean), effective_stds)``""" 
3295          # TODO: how about xcurrent? 
3296          return self.best.get() + ( 
3297              self.countevals, self.countiter, self.gp.pheno(self.mean), self.gp.scales * self.sigma * self.sigma_vec * self.dC**0.5) 
3298   
3299   
3300 -    def result_pretty(self, number_of_runs=0, time_str=None): 
3301          """pretty print result. """ 
3302          s = (' after %i restart' + ('s' if number_of_runs > 1 else '')) \ 
3303              % number_of_runs if number_of_runs else '' 
3304          for k, v in self.stop().items(): 
3305              print('termination on %s=%s%s' % (k, str(v), s) + 
3306                    (' (%s)' % time_str if time_str else '')) 
3307   
3308          print('final/bestever f-value = %e %e' % (self.best.last.f, self.best.f)) 
3309          if self.N < 9: 
3310              print('mean solution: ' + str(self.gp.pheno(self.mean))) 
3311              print('std deviation: ' + str(self.sigma * sqrt(self.dC) * self.gp.scales)) 
3312          else: 
3313              print('mean solution: %s ...]' % (str(self.gp.pheno(self.mean)[:8])[:-1])) 
3314              print('std deviations: %s ...]' % (str((self.sigma * sqrt(self.dC) * self.gp.scales)[:8])[:-1])) 
3315   
3316   
3317 -    def clip_or_fit_solutions(self, pop, idx): 
3318          """make sure that solutions fit to sample distribution, this interface will probably change. 
3319   
3320          In particular the frequency of long vectors appearing in pop[idx] - self.mean is limited. 
3321   
3322          """ 
3323          for k in idx: 
3324              self.repair_genotype(pop[k]) 
3325   
3326 -    def repair_genotype(self, x, copy_if_changed=False): 
3327          """make sure that solutions fit to the sample distribution, this interface will probably change. 
3328   
3329          In particular the frequency of x - self.mean being long is limited. 
3330   
3331          """ 
3332          mold = self.mean 
3333          if 1 < 3:  # hard clip at upper_length 
3334              upper_length = self.N**0.5 + 2 * self.N / (self.N + 2)  # should become an Option, but how? e.g. [0, 2, 2] 
3335              fac = self.mahalanobisNorm(x - mold) / upper_length 
3336   
3337              if fac > 1: 
3338                  if copy_if_changed: 
3339                      x = (x - mold) / fac + mold 
3340                  else:  # should be 25% faster: 
3341                      x -= mold 
3342                      x /= fac 
3343                      x += mold 
3344                  # print self.countiter, k, fac, self.mahalanobisNorm(pop[k] - mold) 
3345                  # adapt also sigma: which are the trust-worthy/injected solutions? 
3346          else: 
3347              if 'checktail' not in self.__dict__:  # hasattr(self, 'checktail') 
3348                  raise NotImplementedError 
3349                  # from check_tail_smooth import CheckTail  # for the time being 
3350                  # self.checktail = CheckTail() 
3351                  # print('untested feature checktail is on') 
3352              fac = self.checktail.addchin(self.mahalanobisNorm(x - mold)) 
3353   
3354              if fac < 1: 
3355                  x = fac * (x - mold) + mold 
3356   
3357          return x 
3358   
3359      # ____________________________________________________________ 
3360      # ____________________________________________________________ 
3361      # 
3362 -    def updateBD(self): 
3363          """update internal variables for sampling the distribution with the 
3364          current covariance matrix C. This method is O(N^3), if C is not diagonal. 
3365   
3366          """ 
3367          # itereigenupdated is always up-to-date in the diagonal case 
3368          # just double check here 
3369          if self.itereigenupdated == self.countiter: 
3370              return 
3371   
3372          # C has already positive updates, here come the additional negative updates 
3373          if self.sp.neg.cmuexp:  # cave: 
3374              if (self.countiter - self.itereigenupdated) * self.sp.neg.cmuexp * self.N < 0.5: # pos.def. guarantied 
3375                  self.C -= self.sp.neg.cmuexp * self.Yneg 
3376              else: # guaranties pos.def. unconditionally 
3377                  # print('exponential update for negative weights (internally more expensive) in iteration', self.countiter) 
3378                  self.update_exponential(self.Yneg, -self.sp.neg.cmuexp) 
3379                  # self.C = self.Ypos + Cs * Mh.expms(-self.sp.neg.cmuexp*Csi*self.Yneg*Csi) * Cs 
3380              self.Yneg = np.zeros((self.N, self.N)) 
3381   
3382          if self.sigma_vec is not 1 and not np.all(self.sigma_vec == 1): 
3383              self.C = dot(dot(np.diag(self.sigma_vec), self.C), np.diag(self.sigma_vec)) 
3384              self.sigma_vec[:] = 1 
3385   
3386          if self.opts['CMA_const_trace'] in (True, 1, 2):  # normalize trace of C 
3387              if self.opts['CMA_const_trace'] == 2: 
3388                  s = np.exp(np.mean(np.log(self.dC))) 
3389              else: 
3390                  s = np.mean(self.dC) 
3391              self.C /= s 
3392              self.dC /= s 
3393          self.C = (self.C + self.C.T) / 2 
3394          # self.C = np.triu(self.C) + np.triu(self.C,1).T  # should work as well 
3395          # self.D, self.B = eigh(self.C) # hermitian, ie symmetric C is assumed 
3396   
3397          if isinstance(self.opts['CMA_eigenmethod'], type(1)): 
3398              print('WARNING: option CMA_eigenmethod should be a function, not an integer') 
3399              if self.opts['CMA_eigenmethod'] == -1: 
3400                  # pygsl 
3401                  # easy to install (well, in Windows install gsl binaries first, 
3402                  # set system path to respective libgsl-0.dll (or cp the dll to 
3403                  # python\DLLS ?), in unzipped pygsl edit 
3404                  # gsl_dist/gsl_site_example.py into gsl_dist/gsl_site.py 
3405                  # and run "python setup.py build" and "python setup.py install" 
3406                  # in MINGW32) 
3407                  if 1 < 3:  # import pygsl on the fly 
3408                      try: 
3409                          import pygsl.eigen.eigenvectors  # TODO efficient enough? 
3410                      except ImportError: 
3411                          print('WARNING: could not find pygsl.eigen module, either install pygsl \n' + 
3412                                '  or set option CMA_eigenmethod=1 (is much slower), option set to 1') 
3413                          self.opts['CMA_eigenmethod'] = 0  # use 0 if 1 is too slow 
3414   
3415                      self.D, self.B = pygsl.eigen.eigenvectors(self.C) 
3416   
3417              elif self.opts['CMA_eigenmethod'] == 0: 
3418                  # TODO: thoroughly test np.linalg.eigh 
3419                  #       numpy.linalg.eig crashes in 200-D 
3420                  #       and EVecs with same EVals are not orthogonal 
3421                  self.D, self.B = np.linalg.eigh(self.C)  # self.B[i] is a row and not an eigenvector 
3422              else:  # is overall two;ten times slower in 10;20-D 
3423                  self.D, self.B = Misc.eig(self.C)  # def eig, see below 
3424          else: 
3425              self.D, self.B = self.opts['CMA_eigenmethod'](self.C) 
3426   
3427   
3428          # assert(sum(self.D-DD) < 1e-6) 
3429          # assert(sum(sum(np.dot(BB, BB.T)-np.eye(self.N))) < 1e-6) 
3430          # assert(sum(sum(np.dot(BB * DD, BB.T) - self.C)) < 1e-6) 
3431          idx = np.argsort(self.D) 
3432          self.D = self.D[idx] 
3433          self.B = self.B[:, idx]  # self.B[i] is a row, columns self.B[:,i] are eigenvectors 
3434          # assert(all(self.B[self.countiter % self.N] == self.B[self.countiter % self.N,:])) 
3435   
3436          # qqqqqqqqqq 
3437          # is O(N^3) 
3438          # assert(sum(abs(self.C - np.dot(self.D * self.B,  self.B.T))) < N**2*1e-11) 
3439   
3440          self.D **= 0.5 
3441          self.itereigenupdated = self.countiter 
3442   
3443 -    def multiplyC(self, alpha): 
3444          """multiply C with a scalar and update all related internal variables (dC, D,...)""" 
3445          self.C *= alpha 
3446          if self.dC is not self.C: 
3447              self.dC *= alpha 
3448          self.D *= alpha**0.5 
3449 -    def update_exponential(self, Z, eta, BDpair=None): 
3450          """exponential update of C that guarantees positive definiteness, that is, 
3451          instead of the assignment ``C = C + eta * Z``, 
3452          we have ``C = C**.5 * exp(eta * C**-.5 * Z * C**-.5) * C**.5``. 
3453   
3454          Parameter `Z` should have expectation zero, e.g. sum(w[i] * z[i] * z[i].T) - C 
3455          if E z z.T = C. 
3456   
3457          Parameter `eta` is the learning rate, for ``eta == 0`` nothing is updated. 
3458   
3459          This function conducts two eigendecompositions, assuming that 
3460          B and D are not up to date, unless `BDpair` is given. Given BDpair, 
3461          B is the eigensystem and D is the vector of sqrt(eigenvalues), one 
3462          eigendecomposition is omitted. 
3463   
3464          Reference: Glasmachers et al 2010, Exponential Natural Evolution Strategies 
3465   
3466          """ 
3467          if eta == 0: 
3468              return 
3469          if BDpair: 
3470              B, D = BDpair 
3471          else: 
3472              D, B = self.opts['CMA_eigenmethod'](self.C) 
3473              D **= 0.5 
3474          Cs = dot(B, (B * D).T)   # square root of C 
3475          Csi = dot(B, (B / D).T)  # square root of inverse of C 
3476          self.C = dot(Cs, dot(Mh.expms(eta * dot(Csi, dot(Z, Csi)), self.opts['CMA_eigenmethod']), Cs)) 
3477   
3478      # ____________________________________________________________ 
3479      # ____________________________________________________________ 
3480 -    def feedForResume(self, X, function_values): 
3481          """Given all "previous" candidate solutions and their respective 
3482          function values, the state of a `CMAEvolutionStrategy` object 
3483          can be reconstructed from this history. This is the purpose of 
3484          function `feedForResume`. 
3485   
3486          Arguments 
3487          --------- 
3488              `X` 
3489                (all) solution points in chronological order, phenotypic 
3490                representation. The number of points must be a multiple 
3491                of popsize. 
3492              `function_values` 
3493                respective objective function values 
3494   
3495          Details 
3496          ------- 
3497          `feedForResume` can be called repeatedly with only parts of 
3498          the history. The part must have the length of a multiple 
3499          of the population size. 
3500          `feedForResume` feeds the history in popsize-chunks into `tell`. 
3501          The state of the random number generator might not be 
3502          reconstructed, but this would be only relevant for the future. 
3503   
3504          Example 
3505          ------- 
3506          :: 
3507   
3508              import cma 
3509   
3510              # prepare 
3511              (x0, sigma0) = ... # initial values from previous trial 
3512              X = ... # list of generated solutions from a previous trial 
3513              f = ... # respective list of f-values 
3514   
3515              # resume 
3516              es = cma.CMAEvolutionStrategy(x0, sigma0) 
3517              es.feedForResume(X, f) 
3518   
3519              # continue with func as objective function 
3520              while not es.stop(): 
3521                 X = es.ask() 
3522                 es.tell(X, [func(x) for x in X]) 
3523   
3524          Credits to Dirk Bueche and Fabrice Marchal for the feeding idea. 
3525   
3526          :See: class `CMAEvolutionStrategy` for a simple dump/load to resume 
3527   
3528          """ 
3529          if self.countiter > 0: 
3530              print('WARNING: feed should generally be used with a new object instance') 
3531          if len(X) != len(function_values): 
3532              raise _Error('number of solutions ' + str(len(X)) + 
3533                  ' and number function values ' + 
3534                  str(len(function_values)) + ' must not differ') 
3535          popsize = self.sp.popsize 
3536          if (len(X) % popsize) != 0: 
3537              raise _Error('number of solutions ' + str(len(X)) + 
3538                      ' must be a multiple of popsize (lambda) ' + 
3539                      str(popsize)) 
3540          for i in rglen((X) / popsize): 
3541              # feed in chunks of size popsize 
3542              self.ask()  # a fake ask, mainly for a conditioned calling of updateBD 
3543                          # and secondary to get possibly the same random state 
3544              self.tell(X[i * popsize:(i + 1) * popsize], function_values[i * popsize:(i + 1) * popsize]) 
3545   
3546      # ____________________________________________________________ 
3547      # ____________________________________________________________ 
3548 -    def readProperties(self): 
3549          """reads dynamic parameters from property file (not implemented) 
3550          """ 
3551          print('not yet implemented') 
3552   
3553      # ____________________________________________________________ 
3554      # ____________________________________________________________ 
3555 -    def mahalanobisNorm(self, dx): 
3556          """compute the Mahalanobis norm that is induced by the adapted sample 
3557          distribution, covariance matrix C times sigma**2. The expected 
3558          Mahalanobis distance to the sample mean is about sqrt(dimension). 
3559   
3560          Argument 
3561          -------- 
3562          A *genotype* difference `dx`. 
3563   
3564          Example 
3565          ------- 
3566          >>> import cma, numpy 
3567          >>> es = cma.CMAEvolutionStrategy(numpy.ones(10), 1) 
3568          >>> xx = numpy.random.randn(2, 10) 
3569          >>> d = es.mahalanobisNorm(es.gp.geno(xx[0]-xx[1])) 
3570   
3571          `d` is the distance "in" the true sample distribution, 
3572          sampled points have a typical distance of ``sqrt(2*es.N)``, 
3573          where `N` is the dimension, and an expected distance of 
3574          close to ``sqrt(N)`` to the sample mean. In the example, 
3575          `d` is the Euclidean distance, because C = I and sigma = 1. 
3576   
3577          """ 
3578          return sqrt(sum((self.D**-1 * np.dot(self.B.T, dx / self.sigma_vec))**2)) / self.sigma 
3579   
3580      # ____________________________________________________________ 
3581      # ____________________________________________________________ 
3582      # 
3583 -    def timesCroot(self, mat): 
3584          """return C**0.5 times mat, where mat can be a vector or matrix. 
3585          Not functional, because _Croot=C**0.5 is never computed (should be in updateBD) 
3586          """ 
3587          print("WARNING: timesCroot is not yet tested") 
3588          if self.opts['CMA_diagonal'] is True \ 
3589                         or self.countiter <= self.opts['CMA_diagonal']: 
3590              res = (self._Croot * mat.T).T 
3591          else: 
3592              res = np.dot(self._Croot, mat) 
3593          return res 
3594 -    def divCroot(self, mat): 
3595          """return C**-1/2 times mat, where mat can be a vector or matrix. 
3596          Not functional, because _Crootinv is never computed. """ 
3597          print("WARNING: divCroot is not yet tested") 
3598          if self.opts['CMA_diagonal'] is True \ 
3599                         or self.countiter <= self.opts['CMA_diagonal']: 
3600              res = (self._Crootinv * mat.T).T 
3601          else: 
3602              res = np.dot(self._Crootinv, mat) 
3603          return res 
3604   
3605      # ____________________________________________________________ 
3606      # ____________________________________________________________ 
3607 -    def disp_annotation(self): 
3608          """print annotation for `disp()`""" 
3609          print('Iterat #Fevals   function value    axis ratio  sigma  minstd maxstd min:sec') 
3610          sys.stdout.flush() 
3611   
3612      # ____________________________________________________________ 
3613      # ____________________________________________________________ 
3614 -    def disp(self, modulo=None):  # TODO: rather assign opt['verb_disp'] as default? 
3615          """prints some infos according to `disp_annotation()`, if 
3616          ``iteration_counter % modulo == 0`` 
3617   
3618          """ 
3619          if modulo is None: 
3620              modulo = self.opts['verb_disp'] 
3621   
3622          # console display 
3623          if modulo: 
3624              if (self.countiter - 1) % (10 * modulo) < 1: 
3625                  self.disp_annotation() 
3626              if self.countiter > 0 and (self.stop() or self.countiter < 4 
3627                                or self.countiter % modulo < 1): 
3628                  if self.opts['verb_time']: 
3629                      toc = self.elapsed_time() 
3630                      stime = str(int(toc // 60)) + ':' + str(round(toc % 60, 1)) 
3631                  else: 
3632                      stime = '' 
3633                  print(' '.join((repr(self.countiter).rjust(5), 
3634                                  repr(self.countevals).rjust(7), 
3635                                  '%.15e' % (min(self.fit.fit)), 
3636                                  '%4.1e' % (self.D.max() / self.D.min()), 
3637                                  '%6.2e' % self.sigma, 
3638                                  '%6.0e' % (self.sigma * sqrt(min(self.dC))), 
3639                                  '%6.0e' % (self.sigma * sqrt(max(self.dC))), 
3640                                  stime))) 
3641                  # if self.countiter < 4: 
3642                  sys.stdout.flush() 
3643   
3644  cma_default_options = { 
3645      # the follow string arguments are evaluated, besides the verb_filenameprefix 
3646      'AdaptSigma': 'CMAAdaptSigmaCSA  # or any other CMAAdaptSigmaBase class e.g. CMAAdaptSigmaTPA', 
3647      'CMA_active': 'True  # negative update, conducted after the original update', 
3648      'CMA_activefac': '1  # learning rate multiplier for active update', 
3649      'CMA_cmean': '1  # learning rate for the mean value', 
3650      'CMA_const_trace': 'False  # normalize trace, value CMA_const_trace=2 normalizes sum log eigenvalues to zero', 
3651      'CMA_diagonal': '0*100*N/sqrt(popsize)  # nb of iterations with diagonal covariance matrix, True for always',  # TODO 4/ccov_separable? 
3652      'CMA_eigenmethod': 'np.linalg.eigh  # 0=numpy-s eigh, -1=pygsl, otherwise cma.Misc.eig (slower)', 
3653      'CMA_elitist': 'False # elitism likely impairs global search performance', 
3654      'CMA_mirrors': 'popsize < 6  # values <0.5 are interpreted as fraction, values >1 as numbers (rounded), otherwise about 0.16 is used', 
3655      'CMA_mu': 'None  # parents selection parameter, default is popsize // 2', 
3656      'CMA_on': 'True  # False or 0 for no adaptation of the covariance matrix', 
3657      'CMA_sample_on_sphere_surface': 'False  #v all mutation vectors have the same length', 
3658      'CMA_rankmu': 'True  # False or 0 for omitting rank-mu update of covariance matrix', 
3659      'CMA_rankmualpha': '0.3  # factor of rank-mu update if mu=1, subject to removal, default might change to 0.0', 
3660      'CMA_dampsvec_fac': 'np.Inf  # tentative and subject to changes, 0.5 would be a "default" damping for sigma vector update', 
3661      'CMA_dampsvec_fade': '0.1  # tentative fading out parameter for sigma vector update', 
3662      'CMA_teststds': 'None  # factors for non-isotropic initial distr. mainly for test purpose, see scaling_...', 
3663      # 'CMA_AII': 'False  # not yet tested', 
3664      'CSA_dampfac': '1  #v positive multiplier for step-size damping, 0.3 is close to optimal on the sphere', 
3665      'CSA_damp_mueff_exponent': '0.5  # zero would mean no dependency of damping on mueff, useful with CSA_disregard_length option', 
3666      'CSA_disregard_length': 'False  #v True is untested', 
3667      'CSA_clip_length_value': 'None  #v untested, [0, 0] means disregarding length completely', 
3668      'CSA_squared': 'False  #v use squared length for sigma-adaptation ', 
3669      'boundary_handling': 'BoundTransform  # or BoundPenalty, unused when ``bounds in (None, [None, None])``', 
3670      'bounds': '[None, None]  # lower (=bounds[0]) and upper domain boundaries, each a scalar or a list/vector', 
3671       # , eval_parallel2': 'not in use {"processes": None, "timeout": 12, "is_feasible": lambda x: True} # distributes function calls to processes processes' 
3672      'fixed_variables': 'None  # dictionary with index-value pairs like {0:1.1, 2:0.1} that are not optimized', 
3673      'ftarget': '-inf  #v target function value, minimization', 
3674      'is_feasible': 'is_feasible  #v a function that computes feasibility, by default lambda x, f: f not in (None, np.NaN)', 
3675      'maxfevals': 'inf  #v maximum number of function evaluations', 
3676      'maxiter': '100 + 50 * (N+3)**2 // popsize**0.5  #v maximum number of iterations', 
3677      'mean_shift_line_samples': 'False #v sample two new solutions colinear to previous mean shift', 
3678      'mindx': '0  #v minimal std in any direction, cave interference with tol*', 
3679      'minstd': '0  #v minimal std in any coordinate direction, cave interference with tol*', 
3680      'pc_line_samples': 'False #v two line samples along the evolution path pc', 
3681      'popsize': '4+int(3*log(N))  # population size, AKA lambda, number of new solution per iteration', 
3682      'randn': 'np.random.standard_normal  #v randn((lam, N)) must return an np.array of shape (lam, N)', 
3683      'scaling_of_variables': 'None  # scale for each variable, sigma0 is interpreted w.r.t. this scale, in that effective_sigma0 = sigma0*scaling. Internally the variables are divided by scaling_of_variables and sigma is unchanged, default is ones(N)', 
3684      'seed': 'None  # random number seed', 
3685      'termination_callback': 'None  #v a function returning True for termination, called after each iteration step and could be abused for side effects', 
3686      'tolfacupx': '1e3  #v termination when step-size increases by tolfacupx (diverges). That is, the initial step-size was chosen far too small and better solutions were found far away from the initial solution x0', 
3687      'tolupsigma': '1e20  #v sigma/sigma0 > tolupsigma * max(sqrt(eivenvals(C))) indicates "creeping behavior" with usually minor improvements', 
3688      'tolfun': '1e-11  #v termination criterion: tolerance in function value, quite useful', 
3689      'tolfunhist': '1e-12  #v termination criterion: tolerance in function value history', 
3690      'tolstagnation': 'int(100 + 100 * N**1.5 / popsize)  #v termination if no improvement over tolstagnation iterations', 
3691      'tolx': '1e-11  #v termination criterion: tolerance in x-changes', 
3692      'transformation': 'None  # [t0, t1] are two mappings, t0 transforms solutions from CMA-representation to f-representation (tf_pheno), t1 is the (optional) back transformation, see class GenoPheno', 
3693      'typical_x': 'None  # used with scaling_of_variables', 
3694      'updatecovwait': 'None  #v number of iterations without distribution update, name is subject to future changes',  # TODO: rename: iterwaitupdatedistribution? 
3695      'verbose': '1  #v verbosity e.v. of initial/final message, -1 is very quiet, not yet fully implemented', 
3696      'verb_append': '0  # initial evaluation counter, if append, do not overwrite output files', 
3697      'verb_disp': '100  #v verbosity: display console output every verb_disp iteration', 
3698      'verb_filenameprefix': 'outcmaes  # output filenames prefix', 
3699      'verb_log': '1  #v verbosity: write data to files every verb_log iteration, writing can be time critical on fast to evaluate functions', 
3700      'verb_plot': '0  #v in fmin(): plot() is called every verb_plot iteration', 
3701      'verb_time': 'True  #v output timings on console', 
3702      'vv': '0  #? versatile variable for hacking purposes, value found in self.opts["vv"]' 
3703  } 
3704 -class CMAOptions(dict): 
3705      """``CMAOptions()`` returns a dictionary with the available options 
3706      and their default values for class ``CMAEvolutionStrategy``. 
3707   
3708      ``CMAOptions('pop')`` returns a subset of recognized options that 
3709      contain 'pop' in there keyword name or (default) value or description. 
3710   
3711      ``CMAOptions(opts)`` returns the subset of recognized options in 
3712      ``dict(opts)``. 
3713   
3714      Option values can be "written" in a string and, when passed to fmin 
3715      or CMAEvolutionStrategy, are evaluated using "N" and "popsize" as 
3716      known values for dimension and population size (sample size, number 
3717      of new solutions per iteration). All default option values are such 
3718      a string. 
3719   
3720      Details 
3721      ------- 
3722      ``CMAOptions`` entries starting with ``tol`` are termination 
3723      "tolerances". 
3724   
3725      For `tolstagnation`, the median over the first and the second half 
3726      of at least `tolstagnation` iterations are compared for both, the 
3727      per-iteration best and per-iteration median function value. 
3728   
3729      Example 
3730      ------- 
3731      :: 
3732   
3733          import cma 
3734          cma.CMAOptions('tol') 
3735   
3736      is a shortcut for cma.CMAOptions().match('tol') that returns all options 
3737      that contain 'tol' in their name or description. 
3738   
3739      :See: `fmin`(), `CMAEvolutionStrategy`, `_CMAParameters` 
3740   
3741      """ 
3742   
3743      # @classmethod # self is the class, not the instance 
3744      # @property 
3745      # def default(self): 
3746      #     """returns all options with defaults""" 
3747      #     return fmin([],[]) 
3748   
3749      @staticmethod 
3750 -    def defaults(): 
3751          """return a dictionary with default option values and description""" 
3752          return dict(cma_default_options) 
3753   
3754      @staticmethod 
3755 -    def versatileOptions(): 
3756          """return list of options that can be changed at any time (not 
3757          only be initialized), however the list might not be entirely up 
3758          to date. 
3759   
3760          The string ' #v ' in the default value indicates a 'versatile' 
3761          option that can be changed any time. 
3762   
3763          """ 
3764          return tuple(sorted(i[0] for i in list(CMAOptions.defaults().items()) if i[1].find(' #v ') > 0)) 
3765   
3766 -    def __init__(self, s=None, unchecked=False): 
3767          """return an `CMAOptions` instance, either with the default 
3768          options, if ``s is None``, or with all options whose name or 
3769          description contains `s`, if `s` is a string (case is 
3770          disregarded), or with entries from dictionary `s` as options, 
3771          not complemented with default options or settings 
3772   
3773          Returns: see above. 
3774   
3775          """ 
3776          # if not CMAOptions.defaults:  # this is different from self.defaults!!! 
3777          #     CMAOptions.defaults = fmin([],[]) 
3778          if s is None: 
3779              super(CMAOptions, self).__init__(CMAOptions.defaults())  # dict.__init__(self, CMAOptions.defaults()) should be the same 
3780              # self = CMAOptions.defaults() 
3781          elif isinstance(s, basestring): 
3782              super(CMAOptions, self).__init__(CMAOptions().match(s)) 
3783              # we could return here 
3784          else: 
3785              super(CMAOptions, self).__init__(s) 
3786   
3787          if not unchecked: 
3788              for key in list(self.keys()): 
3789                  if key not in CMAOptions.defaults(): 
3790                      print('Warning in cma.CMAOptions.__init__(): invalid key ``' + str(key) + '`` removed') 
3791                      self.pop(key) 
3792          # self.evaluated = False  # would become an option entry 
3793   
3794 -    def init(self, dict_or_str, val=None, warn=True): 
3795          """initialize one or several options. 
3796   
3797          Arguments 
3798          --------- 
3799              `dict_or_str` 
3800                  a dictionary if ``val is None``, otherwise a key. 
3801                  If `val` is provided `dict_or_str` must be a valid key. 
3802              `val` 
3803                  value for key 
3804   
3805          Details 
3806          ------- 
3807          Only known keys are accepted. Known keys are in `CMAOptions.defaults()` 
3808   
3809          """ 
3810          # dic = dict_or_key if val is None else {dict_or_key:val} 
3811          dic = dict_or_str 
3812          if val is not None: 
3813              dic = {dict_or_str:val} 
3814   
3815          for key, val in list(dic.items()): 
3816              if key not in CMAOptions.defaults(): 
3817                  # TODO: find a better solution? 
3818                  if warn: 
3819                      print('Warning in cma.CMAOptions.init(): key ' + 
3820                          str(key) + ' ignored') 
3821              else: 
3822                  self[key] = val 
3823   
3824          return self 
3825   
3826 -    def set(self, dic, val=None, warn=True): 
3827          """set can assign versatile options from `CMAOptions.versatileOptions()` 
3828          with a new value, use `init()` for the others. 
3829   
3830          Arguments 
3831          --------- 
3832              `dic` 
3833                  either a dictionary or a key. In the latter 
3834                  case, val must be provided 
3835              `val` 
3836                  value for key 
3837              `warn` 
3838                  bool, print a warning if the option cannot be changed 
3839                  and is therefore omitted 
3840   
3841          This method will be most probably used with the ``opts`` attribute of 
3842          a `CMAEvolutionStrategy` instance. 
3843   
3844          """ 
3845          if val is not None:  # dic is a key in this case 
3846              dic = {dic:val}  # compose a dictionary 
3847          for key, val in list(dic.items()): 
3848              if key in CMAOptions.versatileOptions(): 
3849                  self[key] = val 
3850              elif warn: 
3851                  print('Warning in cma.CMAOptions.set(): key ' + str(key) + ' ignored') 
3852          return self  # to allow o = CMAOptions(o).set(new) 
3853   
3854 -    def complement(self): 
3855          """add all missing options with their default values""" 
3856   
3857          for key in CMAOptions.defaults(): 
3858              if key not in self: 
3859                  self[key] = CMAOptions.defaults()[key] 
3860          return self 
3861   
3862 -    def settable(self): 
3863          """return the subset of those options that are settable at any 
3864          time. 
3865   
3866          Settable options are in `versatileOptions()`, but the 
3867          list might be incomplete. 
3868   
3869          """ 
3870          return CMAOptions([i for i in list(self.items()) 
3871                                  if i[0] in CMAOptions.versatileOptions()]) 
3872   
3873 -    def __call__(self, key, default=None, loc=None): 
3874          """evaluate and return the value of option `key` on the fly, or 
3875          returns those options whose name or description contains `key`, 
3876          case disregarded. 
3877   
3878          Details 
3879          ------- 
3880          Keys that contain `filename` are not evaluated. 
3881          For ``loc==None``, `self` is used as environment 
3882          but this does not define `N`. 
3883   
3884          :See: `eval()`, `evalall()` 
3885   
3886          """ 
3887          try: 
3888              val = self[key] 
3889          except: 
3890              return self.match(key) 
3891   
3892          if loc is None: 
3893              loc = self  # TODO: this hack is not so useful: popsize could be there, but N is missing 
3894          try: 
3895              if isinstance(val, basestring): 
3896                  val = val.split('#')[0].strip()  # remove comments 
3897                  if isinstance(val, basestring) and key.find('filename') < 0 and key.find('mindx') < 0: 
3898                      val = eval(val, globals(), loc) 
3899              # invoke default 
3900              # TODO: val in ... fails with array type, because it is applied element wise! 
3901              # elif val in (None,(),[],{}) and default is not None: 
3902              elif val is None and default is not None: 
3903                  val = eval(str(default), globals(), loc) 
3904          except: 
3905              pass  # slighly optimistic: the previous is bug-free 
3906          return val 
3907   
3908 -    def eval(self, key, default=None, loc=None): 
3909          """Evaluates and sets the specified option value in 
3910          environment `loc`. Many options need `N` to be defined in 
3911          `loc`, some need `popsize`. 
3912   
3913          Details 
3914          ------- 
3915          Keys that contain 'filename' are not evaluated. 
3916          For `loc` is None, the self-dict is used as environment 
3917   
3918          :See: `evalall()`, `__call__` 
3919   
3920          """ 
3921          # TODO: try: loc['dim'] = loc['N'] etc 
3922          self[key] = self(key, default, loc) 
3923          return self[key] 
3924   
3925 -    def evalall(self, loc=None, defaults=None): 
3926          """Evaluates all option values in environment `loc`. 
3927   
3928          :See: `eval()` 
3929   
3930          """ 
3931          if defaults is None: 
3932              defaults = cma_default_options 
3933          # TODO: this needs rather the parameter N instead of loc 
3934          if 'N' in loc:  # TODO: __init__ of CMA can be simplified 
3935              popsize = self('popsize', defaults['popsize'], loc) 
3936              for k in list(self.keys()): 
3937                  self.eval(k, defaults[k], 
3938                            {'N':loc['N'], 'popsize':popsize}) 
3939          return self 
3940   
3941 -    def match(self, s=''): 
3942          """return all options that match, in the name or the description, 
3943          with string `s`, case is disregarded. 
3944   
3945          Example: ``cma.CMAOptions().match('verb')`` returns the verbosity options. 
3946   
3947          """ 
3948          match = s.lower() 
3949          res = {} 
3950          for k in sorted(self): 
3951              s = str(k) + '=\'' + str(self[k]) + '\'' 
3952              if match in s.lower(): 
3953                  res[k] = self[k] 
3954          return CMAOptions(res) 
3955   
3956 -    def pp(self): 
3957          pprint(self) 
3958   
3959 -    def pprint(self, linebreak=80): 
3960          for i in sorted(self.items()): 
3961              s = str(i[0]) + "='" + str(i[1]) + "'" 
3962              a = s.split(' ') 
3963   
3964              # print s in chunks 
3965              l = ''  # start entire to the left 
3966              while a: 
3967                  while a and len(l) + len(a[0]) < linebreak: 
3968                      l += ' ' + a.pop(0) 
3969                  print(l) 
3970                  l = '        '  # tab for subsequent lines 
3971      print_ = pprint  # Python style to prevent clash with keywords 
3972      printme = pprint 
3973   
3974  # ____________________________________________________________ 
3975  # ____________________________________________________________ 
3976 -class CMAStopDict(dict): 
3977      """keep and update a termination condition dictionary, which is 
3978      "usually" empty and returned by `CMAEvolutionStrategy.stop()`. 
3979      The class methods entirely depend on `CMAEvolutionStrategy` class 
3980      attributes. 
3981   
3982      Details 
3983      ------- 
3984      This class is not relevant for the end-user and could be a nested 
3985      class, but nested classes cannot be serialized. 
3986   
3987      Example 
3988      ------- 
3989      >>> import cma 
3990      >>> sd = cma.CMAStopDict() 
3991      >>> es = cma.CMAEvolutionStrategy(4 * [1], 1, {'verbose':-1}) 
3992      >>> print(sd(es)) 
3993      {} 
3994      >>> es.optimize(cma.fcts.sphere, verb_disp=0) 
3995      >>> print(sd(es)) 
3996      {'tolfun': 1e-11} 
3997      >>> assert sd(es) == es.stop() 
3998   
3999      :See: `OOOptimizer.stop()`, `CMAEvolutionStrategy.stop()` 
4000   
4001      """ 
4002 -    def __init__(self, d={}): 
4003          update = isinstance(d, CMAEvolutionStrategy) 
4004          inherit = isinstance(d, CMAStopDict) 
4005          super(CMAStopDict, self).__init__({} if update else d) 
4006          self._stoplist = d._stoplist if inherit else []  # multiple entries 
4007          self.lastiter = d.lastiter if inherit else 0  # probably not necessary 
4008          if update: 
4009              self._update(d) 
4010   
4011 -    def __call__(self, es=None): 
4012          """update and return the termination conditions dictionary 
4013   
4014          """ 
4015          if es is None and self.es is None: 
4016              raise ValueError('termination conditions need an optimizer to act upon') 
4017          self._update(es) 
4018          return self 
4019   
4020 -    def _update(self, es): 
4021          """Test termination criteria and update dictionary 
4022   
4023          """ 
4024          if es is None: 
4025              es = self.es 
4026          assert es is not None 
4027          if es.countiter == self.lastiter: 
4028              if es.countiter == 0: 
4029                  self.__init__() 
4030                  return self 
4031              try: 
4032                  if es == self.es: 
4033                      return self 
4034              except:  # self.es not yet assigned 
4035                  pass 
4036   
4037          self.lastiter = es.countiter 
4038          self.es = es 
4039   
4040          self.stoplist = [] 
4041   
4042          N = es.N 
4043          opts = es.opts 
4044          self.opts = opts  # a hack to get _addstop going 
4045   
4046          # fitness: generic criterion, user defined w/o default 
4047          self._addstop('ftarget', 
4048                        es.best.f < opts['ftarget']) 
4049          # maxiter, maxfevals: generic criteria 
4050          self._addstop('maxfevals', 
4051                        es.countevals - 1 >= opts['maxfevals']) 
4052          self._addstop('maxiter', 
4053                        es.countiter >= opts['maxiter']) 
4054          # tolx, tolfacupx: generic criteria 
4055          # tolfun, tolfunhist (CEC:tolfun includes hist) 
4056          self._addstop('tolx', 
4057                        all([es.sigma * xi < opts['tolx'] for xi in es.pc]) and 
4058                        all([es.sigma * xi < opts['tolx'] for xi in sqrt(es.dC)])) 
4059          self._addstop('tolfacupx', 
4060                        any([es.sigma * sig > es.sigma0 * opts['tolfacupx'] 
4061                             for sig in sqrt(es.dC)])) 
4062          self._addstop('tolfun', 
4063                        es.fit.fit[-1] - es.fit.fit[0] < opts['tolfun'] and 
4064                        max(es.fit.hist) - min(es.fit.hist) < opts['tolfun']) 
4065          self._addstop('tolfunhist', 
4066                        len(es.fit.hist) > 9 and 
4067                        max(es.fit.hist) - min(es.fit.hist) < opts['tolfunhist']) 
4068   
4069          # worst seen false positive: table N=80,lam=80, getting worse for fevals=35e3 \approx 50 * N**1.5 
4070          # but the median is not so much getting worse 
4071          # / 5 reflects the sparsity of histbest/median 
4072          # / 2 reflects the left and right part to be compared 
4073          l = int(max((opts['tolstagnation'] / 5. / 2, len(es.fit.histbest) / 10))) 
4074          # TODO: why max(..., len(histbest)/10) ??? 
4075          # TODO: the problem in the beginning is only with best ==> ??? 
4076          # equality should handle flat fitness 
4077          self._addstop('tolstagnation',  # leads sometimes early stop on ftablet, fcigtab, N>=50? 
4078                        1 < 3 and opts['tolstagnation'] and es.countiter > N * (5 + 100 / es.popsize) and 
4079                        len(es.fit.histbest) > 100 and 2 * l < len(es.fit.histbest) and 
4080                        np.median(es.fit.histmedian[:l]) >= np.median(es.fit.histmedian[l:2 * l]) and 
4081                        np.median(es.fit.histbest[:l]) >= np.median(es.fit.histbest[l:2 * l])) 
4082          # iiinteger: stagnation termination can prevent to find the optimum 
4083   
4084          self._addstop('tolupsigma', opts['tolupsigma'] and 
4085                        es.sigma / es.sigma0 / np.max(es.D) > opts['tolupsigma']) 
4086   
4087          if 1 < 3: 
4088              # non-user defined, method specific 
4089              # noeffectaxis (CEC: 0.1sigma), noeffectcoord (CEC:0.2sigma), conditioncov 
4090              self._addstop('noeffectcoord', 
4091                           any([es.mean[i] == es.mean[i] + 0.2 * es.sigma * sqrt(es.dC[i]) 
4092                                for i in xrange(N)])) 
4093              if opts['CMA_diagonal'] is not True and es.countiter > opts['CMA_diagonal']: 
4094                  i = es.countiter % N 
4095                  self._addstop('noeffectaxis', 
4096                               sum(es.mean == es.mean + 0.1 * es.sigma * es.D[i] * es.B[:, i]) == N) 
4097              self._addstop('conditioncov', 
4098                           es.D[-1] > 1e7 * es.D[0], 1e14)  # TODO 
4099   
4100              self._addstop('callback', es.callbackstop)  # termination_callback 
4101          if len(self): 
4102              self._addstop('flat fitness: please (re)consider how to compute the fitness more elaborate', 
4103                           len(es.fit.hist) > 9 and 
4104                           max(es.fit.hist) == min(es.fit.hist)) 
4105          return self 
4106   
4107 -    def _addstop(self, key, cond, val=None): 
4108          if cond: 
4109              self.stoplist.append(key)  # can have the same key twice 
4110              self[key] = self.opts.get(key, None) 
4111   
4112 -    def clear(self): 
4113          for k in list(self): 
4114              self.pop(k) 
4115          self.stoplist = [] 
4116   
4117  # ____________________________________________________________ 
4118  # ____________________________________________________________ 
4119 -class _CMAParameters(object): 
4120      """strategy parameters like population size and learning rates. 
4121   
4122      Note: 
4123          contrary to `CMAOptions`, `_CMAParameters` is not (yet) part of the 
4124          "user-interface" and subject to future changes (it might become 
4125          a `collections.namedtuple`) 
4126   
4127      Example 
4128      ------- 
4129      >>> import cma 
4130      >>> es = cma.CMAEvolutionStrategy(20 * [0.1], 1) 
4131      (6_w,12)-CMA-ES (mu_w=3.7,w_1=40%) in dimension 20 (seed=504519190)  # the seed is "random" by default 
4132      >>> 
4133      >>> type(es.sp)  # sp contains the strategy parameters 
4134      <class 'cma._CMAParameters'> 
4135      >>> 
4136      >>> es.sp.disp() 
4137      {'CMA_on': True, 
4138       'N': 20, 
4139       'c1': 0.004181139918745593, 
4140       'c1_sep': 0.034327992810300939, 
4141       'cc': 0.17176721127681213, 
4142       'cc_sep': 0.25259494835857677, 
4143       'cmean': 1.0, 
4144       'cmu': 0.0085149624979034746, 
4145       'cmu_sep': 0.057796356229390715, 
4146       'cs': 0.21434997799189287, 
4147       'damps': 1.2143499779918929, 
4148       'mu': 6, 
4149       'mu_f': 6.0, 
4150       'mueff': 3.7294589343030671, 
4151       'popsize': 12, 
4152       'rankmualpha': 0.3, 
4153       'weights': array([ 0.40240294,  0.25338908,  0.16622156,  0.10437523,  0.05640348, 
4154              0.01720771])} 
4155      >>> 
4156      >> es.sp == cma._CMAParameters(20, 12, cma.CMAOptions().evalall({'N': 20})) 
4157      True 
4158   
4159      :See: `CMAOptions`, `CMAEvolutionStrategy` 
4160   
4161      """ 
4162 -    def __init__(self, N, opts, ccovfac=1, verbose=True): 
4163          """Compute strategy parameters, mainly depending on 
4164          dimension and population size, by calling `set` 
4165   
4166          """ 
4167          self.N = N 
4168          if ccovfac == 1: 
4169              ccovfac = opts['CMA_on']  # that's a hack 
4170          self.popsize = None  # declaring the attribute, not necessary though 
4171          self.set(opts, ccovfac=ccovfac, verbose=verbose) 
4172   
4173 -    def set(self, opts, popsize=None, ccovfac=1, verbose=True): 
4174          """Compute strategy parameters as a function 
4175          of dimension and population size """ 
4176   
4177          alpha_cc = 1.0  # cc-correction for mueff, was zero before 
4178   
4179          def cone(df, mu, N, alphacov=2.0): 
4180              """rank one update learning rate, ``df`` is disregarded and obsolete, reduce alphacov on noisy problems, say to 0.5""" 
4181              return alphacov / ((N + 1.3)**2 + mu) 
4182   
4183          def cmu(df, mu, alphamu=0.0, alphacov=2.0): 
4184              """rank mu learning rate, disregarding the constrant cmu <= 1 - cone""" 
4185              c = alphacov * (alphamu + mu - 2 + 1 / mu) / ((N + 2)**2 + alphacov * mu / 2) 
4186              # c = alphacov * (alphamu + mu - 2 + 1/mu) / (2 * (N + 2)**1.5 + alphacov * mu / 2) 
4187              # print 'cmu =', c 
4188              return c 
4189   
4190          def conedf(df, mu, N): 
4191              """used for computing separable learning rate""" 
4192              return 1. / (df + 2.*sqrt(df) + float(mu) / N) 
4193   
4194          def cmudf(df, mu, alphamu): 
4195              """used for computing separable learning rate""" 
4196              return (alphamu + mu - 2. + 1. / mu) / (df + 4.*sqrt(df) + mu / 2.) 
4197   
4198          sp = self 
4199          N = sp.N 
4200          if popsize: 
4201              opts.evalall({'N':N, 'popsize':popsize}) 
4202          else: 
4203              popsize = opts.evalall({'N':N})['popsize']  # the default popsize is computed in CMAOptions() 
4204          sp.popsize = popsize 
4205          if opts['CMA_mirrors'] < 0.5: 
4206              sp.lam_mirr = int(0.5 + opts['CMA_mirrors'] * popsize) 
4207          elif opts['CMA_mirrors'] > 1: 
4208              sp.lam_mirr = int(0.5 + opts['CMA_mirrors']) 
4209          else: 
4210              sp.lam_mirr = int(0.5 + 0.16 * min((popsize, 2 * N + 2)) + 0.29)  # 0.158650... * popsize is optimal 
4211              # lam = arange(2,22) 
4212              # mirr = 0.16 + 0.29/lam 
4213              # print(lam); print([int(0.5 + l) for l in mirr*lam]) 
4214              # [ 2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21] 
4215              # [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4] 
4216   
4217          sp.mu_f = sp.popsize / 2.0  # float value of mu 
4218          if opts['CMA_mu'] is not None: 
4219              sp.mu_f = opts['CMA_mu'] 
4220          sp.mu = int(sp.mu_f + 0.499999)  # round down for x.5 
4221          # in principle we have mu_opt = popsize/2 + lam_mirr/2, 
4222          # which means in particular weights should only be negative for q > 0.5+mirr_frac/2 
4223          if sp.mu > sp.popsize - 2 * sp.lam_mirr + 1: 
4224              print("WARNING: pairwise selection is not implemented, therefore " + 
4225                    " mu = %d > %d = %d - 2*%d + 1 = popsize - 2*mirr + 1 can produce a bias" % ( 
4226                      sp.mu, sp.popsize - 2 * sp.lam_mirr + 1, sp.popsize, sp.lam_mirr)) 
4227          if sp.lam_mirr > sp.popsize // 2: 
4228              raise _Error("fraction of mirrors in the population as read from option CMA_mirrors cannot be larger 0.5, " + 
4229                           "theoretically optimal is 0.159") 
4230          sp.weights = log(max([sp.mu, sp.popsize / 2.0]) + 0.5) - log(1 + np.arange(sp.mu)) 
4231          sp.weights /= sum(sp.weights) 
4232          sp.mueff = 1 / sum(sp.weights**2) 
4233          sp.cs = (sp.mueff + 2) / (N + sp.mueff + 3) 
4234          # TODO: clean up (here the cumulation constant is shorter if sigma_vec is used) 
4235          sp.dampsvec = opts['CMA_dampsvec_fac'] * (N + 2) if opts['CMA_dampsvec_fac'] else np.Inf 
4236          sp.dampsvec_fading = opts['CMA_dampsvec_fade'] 
4237          if np.isfinite(sp.dampsvec): 
4238              sp.cs = ((sp.mueff + 2) / (N + sp.mueff + 3))**0.5 
4239          # sp.cs = (sp.mueff + 2) / (N + 1.5*sp.mueff + 1) 
4240          sp.cc = (4 + alpha_cc * sp.mueff / N) / (N + 4 + alpha_cc * 2 * sp.mueff / N) 
4241          sp.cc_sep = (1 + 1 / N + alpha_cc * sp.mueff / N) / (N**0.5 + 1 / N + alpha_cc * 2 * sp.mueff / N)  # \not\gg\cc 
4242          sp.rankmualpha = opts['CMA_rankmualpha'] 
4243          # sp.rankmualpha = _evalOption(opts['CMA_rankmualpha'], 0.3) 
4244          sp.c1 = ccovfac * min(1, sp.popsize / 6) * cone((N**2 + N) / 2, sp.mueff, N)  # 2. / ((N+1.3)**2 + sp.mucov) 
4245          sp.c1_sep = ccovfac * conedf(N, sp.mueff, N) 
4246          if opts['CMA_rankmu'] != 0:  # also empty 
4247              sp.cmu = min(1 - sp.c1, ccovfac * cmu((N**2 + N) / 2, sp.mueff, sp.rankmualpha)) 
4248              sp.cmu_sep = min(1 - sp.c1_sep, ccovfac * cmudf(N, sp.mueff, sp.rankmualpha)) 
4249          else: 
4250              sp.cmu = sp.cmu_sep = 0 
4251   
4252          sp.neg = _BlancClass() 
4253          if opts['CMA_active']: 
4254              # in principle we have mu_opt = popsize/2 + lam_mirr/2, 
4255              # which means in particular weights should only be negative for q > 0.5+mirr_frac/2 
4256              if 1 < 3: # seems most natural: continuation of log(lambda/2) - log(k) qqqqqqqqqqqqqqqqqqqqqqqqqq 
4257                  sp.neg.mu_f = popsize // 2  # not sure anymore what this is good for 
4258                  sp.neg.weights = array([log(k) - log(popsize/2 + 1/2) for k in np.arange(np.ceil(popsize/2 + 1.1/2), popsize + .1)]) 
4259              sp.neg.mu = len(sp.neg.weights) 
4260              sp.neg.weights /= sum(sp.neg.weights) 
4261              sp.neg.mueff = 1 / sum(sp.neg.weights**2) 
4262              sp.neg.cmuexp = opts['CMA_activefac'] * 0.5 * sp.neg.mueff / ((N + 2)**1.5 + 1.0 * sp.neg.mueff) 
4263              # reasoning on learning rate cmuexp: with sum |w| == 1 and 
4264              #   length-normalized vectors in the update, the residual 
4265              #   variance in any direction exceeds exp(-N*cmuexp) 
4266              assert sp.neg.mu >= sp.lam_mirr  # not really necessary 
4267              # sp.neg.minresidualvariance = 0.66  # not it use, keep at least 0.66 in all directions, small popsize is most critical 
4268          else: 
4269              sp.neg.cmuexp = 0 
4270   
4271          sp.CMA_on = sp.c1 + sp.cmu > 0 
4272          # print(sp.c1_sep / sp.cc_sep) 
4273   
4274          if not opts['CMA_on'] and opts['CMA_on'] not in (None, [], (), ''): 
4275              sp.CMA_on = False 
4276              # sp.c1 = sp.cmu = sp.c1_sep = sp.cmu_sep = 0 
4277          mueff_exponent = 0.5 
4278          if 1 < 3: 
4279              mueff_exponent = opts['CSA_damp_mueff_exponent'] 
4280          # TODO: this will disappear, as it is done in class CMAAdaptSigmaCSA 
4281          sp.damps = opts['CSA_dampfac'] * (0.5 + 
4282                                            0.5 * min([1, (sp.lam_mirr / (0.159 * sp.popsize) - 1)**2])**1 + 
4283                                            2 * max([0, ((sp.mueff - 1) / (N + 1))**mueff_exponent - 1]) + sp.cs 
4284                                            ) 
4285          sp.cmean = float(opts['CMA_cmean']) 
4286          # sp.kappa = 1  # 4-D, lam=16, rank1, kappa < 4 does not influence convergence rate 
4287                          # in larger dim it does, 15-D with defaults, kappa=8 factor 2 
4288          if sp.cmean != 1: 
4289              print('  cmean = %f' % (sp.cmean)) 
4290   
4291          if verbose: 
4292              if not sp.CMA_on: 
4293                  print('covariance matrix adaptation turned off') 
4294              if opts['CMA_mu'] != None: 
4295                  print('mu = %f' % (sp.mu_f)) 
4296   
4297          # return self  # the constructor returns itself 
4298   
4299 -    def disp(self): 
4300          pprint(self.__dict__) 
4301   
4302 -def fmin(objective_function, x0, sigma0, 
4303           options=None, 
4304           args=(), 
4305           restarts=0, 
4306           restart_from_best='False', 
4307           incpopsize=2, 
4308           eval_initial_x=False, 
4309           noise_handler=None, 
4310           noise_change_sigma_exponent=1, 
4311           noise_kappa_exponent=0  # TODO: add max kappa value as parameter 
4312          ): 
4313      """functional interface to the stochastic optimizer CMA-ES 
4314      for non-convex function minimization. 
4315   
4316      Calling Sequences 
4317      ================= 
4318          ``fmin(objective_function, x0, sigma0)`` 
4319              minimizes `objective_function` starting at `x0` and with standard deviation 
4320              `sigma0` (step-size) 
4321          ``fmin(objective_function, x0, sigma0, options={'ftarget': 1e-5})`` 
4322              minimizes `objective_function` up to target function value 1e-5, which 
4323              is typically useful for benchmarking. 
4324          ``fmin(objective_function, x0, sigma0, args=('f',))`` 
4325              minimizes `objective_function` called with an additional argument ``'f'``. 
4326          ``fmin(objective_function, x0, sigma0, options={'ftarget':1e-5, 'popsize':40})`` 
4327              uses additional options ``ftarget`` and ``popsize`` 
4328          ``fmin(objective_function, esobj, None, options={'maxfevals': 1e5})`` 
4329              uses the `CMAEvolutionStrategy` object instance `esobj` to optimize 
4330              `objective_function`, similar to `esobj.optimize()`. 
4331   
4332      Arguments 
4333      ========= 
4334          `objective_function` 
4335              function to be minimized. Called as ``objective_function(x,*args)``. 
4336              `x` is a one-dimensional `numpy.ndarray`. `objective_function` 
4337              can return `numpy.NaN`, 
4338              which is interpreted as outright rejection of solution `x` 
4339              and invokes an immediate resampling and (re-)evaluation 
4340              of a new solution not counting as function evaluation. 
4341          `x0` 
4342              list or `numpy.ndarray`, initial guess of minimum solution 
4343              before the application of the geno-phenotype transformation 
4344              according to the ``transformation`` option. Otherwise 
4345              `x0` can also be a `cma.CMAEvolutionStrategy` object instance. 
4346              In the latter case `sigma0` can be ``None``. 
4347          `sigma0` 
4348              scalar, initial standard deviation in each coordinate. 
4349              `sigma0` should be about 1/4th of the search domain width (where the 
4350              optimum is to be expected). The variables in `objective_function` 
4351              should be scaled such that they presumably have similar sensitivity. 
4352              See also option `scaling_of_variables`. 
4353          `options` 
4354              a dictionary with additional options passed to the constructor 
4355              of class ``CMAEvolutionStrategy``, see ``cma.CMAOptions()`` for 
4356              a list of available options. 
4357          ``args=()`` 
4358              arguments to be used to call the `objective_function` 
4359          ``restarts=0`` 
4360              number of restarts 
4361          ``restart_from_best=False`` 
4362              which point to restart from 
4363          ``incpopsize=2`` 
4364              multiplier for increasing the population size `popsize` before each restart 
4365          ``eval_initial_x=False`` 
4366              evaluate initial solution 
4367          ``noise_handler=None`` 
4368              a ``class NoiseHandler`` object or ``None`` 
4369          ``noise_change_sigma_exponent=1`` 
4370              exponent for sigma increment for additional noise treatment 
4371          ``noise_evaluations_as_kappa`` 
4372              instead of applying reevaluations, the "number of evaluations" 
4373              is (ab)used as scaling factor kappa (experimental). 
4374   
4375      Optional Arguments 
4376      ================== 
4377      All values in the `options` dictionary are evaluated if they are of 
4378      type `str`, besides `verb_filenameprefix`, see class `CMAOptions` for details. 
4379      The full list is available in ``cma.default_options``. 
4380   
4381      >>> import cma 
4382      >>> cma.CMAOptions() 
4383   
4384      Subsets of options can be displayed, for example like ``cma.CMAOptions('tol')``, 
4385      or ``cma.CMAOptions('bound')``, see also class `CMAOptions`. 
4386   
4387      Return 
4388      ====== 
4389      Similar to `OOOptimizer.optimize()` and/or `CMAEvolutionStrategy.optimize()`, return the 
4390      list provided by `CMAEvolutionStrategy.result()` appended with an `OOOptimizer` and an 
4391      `BaseDataLogger`:: 
4392   
4393          res = es.result() + (es.stop(), es, logger) 
4394   
4395      where 
4396          - ``res[0]`` (``xopt``) -- best evaluated solution 
4397          - ``res[1]`` (``fopt``) -- respective function value 
4398          - ``res[2]`` (``evalsopt``) -- respective number of function evaluations 
4399          - ``res[3]`` (``evals``) -- number of overall conducted objective function evaluations 
4400          - ``res[4]`` (``iterations``) -- number of overall conducted iterations 
4401          - ``res[5]`` (``xmean``) -- mean of the final sample distribution 
4402          - ``res[6]`` (``stds``) -- effective stds of the final sample distribution 
4403          - ``res[-3]`` (``stop``) -- termination condition(s) in a dictionary 
4404          - ``res[-2]`` (``cmaes``) -- class `CMAEvolutionStrategy` instance 
4405          - ``res[-1]`` (``logger``) -- class `CMADataLogger` instance 
4406   
4407      Details 
4408      ======= 
4409      This function is an interface to the class `CMAEvolutionStrategy`. The 
4410      latter class should be used when full control over the iteration loop 
4411      of the optimizer is desired. 
4412   
4413      The noise handling follows closely [Hansen et al 2009, A Method for Handling 
4414      Uncertainty in Evolutionary Optimization...] in the measurement part, but the 
4415      implemented treatment is slightly different: for ``noiseS > 0``, ``evaluations`` 
4416      (time) and sigma are increased by ``alpha``. For ``noiseS < 0``, ``evaluations`` 
4417      (time) is decreased by ``alpha**(1/4)``. The option ``noise_handling`` switches 
4418      the noise handling on/off, the given value defines the maximal number 
4419      of evaluations for a single fitness computation. If ``noise_handling`` is a list, 
4420      the smallest element defines the minimal number and if the list has three elements, 
4421      the median value is the start value for ``evaluations``. See also class 
4422      `NoiseHandler`. 
4423   
4424      Examples 
4425      ======== 
4426      The following example calls `fmin` optimizing the Rosenbrock function 
4427      in 10-D with initial solution 0.1 and initial step-size 0.5. The 
4428      options are specified for the usage with the `doctest` module. 
4429   
4430      >>> import cma 
4431      >>> # cma.CMAOptions()  # returns all possible options 
4432      >>> options = {'CMA_diagonal':100, 'seed':1234, 'verb_time':0} 
4433      >>> 
4434      >>> res = cma.fmin(cma.fcts.rosen, [0.1] * 10, 0.5, options) 
4435      (5_w,10)-CMA-ES (mu_w=3.2,w_1=45%) in dimension 10 (seed=1234) 
4436         Covariance matrix is diagonal for 10 iterations (1/ccov=29.0) 
4437      Iterat #Fevals   function value     axis ratio  sigma   minstd maxstd min:sec 
4438          1      10 1.264232686260072e+02 1.1e+00 4.40e-01  4e-01  4e-01 
4439          2      20 1.023929748193649e+02 1.1e+00 4.00e-01  4e-01  4e-01 
4440          3      30 1.214724267489674e+02 1.2e+00 3.70e-01  3e-01  4e-01 
4441        100    1000 6.366683525319511e+00 6.2e+00 2.49e-02  9e-03  3e-02 
4442        200    2000 3.347312410388666e+00 1.2e+01 4.52e-02  8e-03  4e-02 
4443        300    3000 1.027509686232270e+00 1.3e+01 2.85e-02  5e-03  2e-02 
4444        400    4000 1.279649321170636e-01 2.3e+01 3.53e-02  3e-03  3e-02 
4445        500    5000 4.302636076186532e-04 4.6e+01 4.78e-03  3e-04  5e-03 
4446        600    6000 6.943669235595049e-11 5.1e+01 5.41e-06  1e-07  4e-06 
4447        650    6500 5.557961334063003e-14 5.4e+01 1.88e-07  4e-09  1e-07 
4448      termination on tolfun : 1e-11 
4449      final/bestever f-value = 5.55796133406e-14 2.62435631419e-14 
4450      mean solution:  [ 1.          1.00000001  1.          1. 
4451          1.          1.00000001  1.00000002  1.00000003 ...] 
4452      std deviation: [ 3.9193387e-09  3.7792732e-09  4.0062285e-09  4.6605925e-09 
4453          5.4966188e-09   7.4377745e-09   1.3797207e-08   2.6020765e-08 ...] 
4454      >>> 
4455      >>> print('best solutions fitness = %f' % (res[1])) 
4456      best solutions fitness = 2.62435631419e-14 
4457      >>> assert res[1] < 1e-12 
4458   
4459      The above call is pretty much equivalent with the slightly more 
4460      verbose call :: 
4461   
4462          res = cma.CMAEvolutionStrategy([0.1] * 10, 0.5, 
4463                      options=options).optimize(cma.fcts.rosen) 
4464   
4465      In either case, the method :: 
4466   
4467          cma.plot(); 
4468   
4469      (based on `matplotlib.pyplot`) produces a plot of the run and, if 
4470      necessary:: 
4471   
4472          cma.show() 
4473   
4474      shows the plot in a window. To continue you might need to 
4475      close the pop-up window. This behavior seems to disappear in 
4476      subsequent calls of `cma.plot()` and is avoided by using 
4477      `ipython` with `-pylab` option. Finally :: 
4478   
4479          cma.savefig('myfirstrun')  # savefig from matplotlib.pyplot 
4480   
4481      will save the figure in a png. 
4482   
4483      :See: `CMAEvolutionStrategy`, `OOOptimizer.optimize(), `plot()`, 
4484          `CMAOptions`, `scipy.optimize.fmin()` 
4485   
4486      """  # style guides say there should be the above empty line 
4487      if 1 < 3:  # try: # pass on KeyboardInterrupt 
4488          if not objective_function:  # return available options in a dictionary 
4489              return CMAOptions()  # these opts are by definition valid 
4490   
4491          fmin_options = locals().copy()  # archive original options 
4492          del fmin_options['objective_function'] 
4493          del fmin_options['x0'] 
4494          del fmin_options['sigma0'] 
4495          del fmin_options['options'] 
4496          del fmin_options['args'] 
4497   
4498          if options is None: 
4499              options = cma_default_options 
4500          opts = CMAOptions(options.copy()).complement() 
4501   
4502          irun = 0 
4503          best = BestSolution() 
4504          while True:  # restart loop 
4505              # recover from a CMA object 
4506              if irun == 0 and isinstance(x0, CMAEvolutionStrategy): 
4507                  es = x0 
4508                  x0 = es.inputargs['x0']  # for the next restarts 
4509                  if sigma0 is None or not np.isscalar(array(sigma0)): 
4510                      sigma0 = es.inputargs['sigma0']  # for the next restarts 
4511                  # ignore further input args and keep original options 
4512              else:  # default case 
4513                  if irun and eval(str(fmin_options['restart_from_best'])): 
4514                      print('CAVE: restart_from_best is often not useful') 
4515                      es = CMAEvolutionStrategy(best.x, sigma0, opts) 
4516                  else: 
4517                      es = CMAEvolutionStrategy(x0, sigma0, opts) 
4518                  if eval_initial_x: 
4519                      x = es.gp.pheno(es.mean, into_bounds=es.boundary_handler.repair, archive=es.sent_solutions) 
4520                      es.best.update([x], es.sent_solutions, [objective_function(x, *args)], 1) 
4521                      es.countevals += 1 
4522   
4523              opts = es.opts  # processed options, unambiguous 
4524              # a hack: 
4525              fmin_opts = CMAOptions(fmin_options.copy(), unchecked=True) 
4526              for k in fmin_opts: 
4527                  # locals() cannot be modified directly, exec won't work in 3.x, therefore 
4528                  fmin_opts.eval(k, loc={'N': es.N, 'popsize': opts['popsize']}) 
4529   
4530              append = opts['verb_append'] or es.countiter > 0 or irun > 0 
4531              # es.logger is "the same" logger, because the "identity" is only determined by the `filenameprefix` 
4532              logger = CMADataLogger(opts['verb_filenameprefix'], opts['verb_log']) 
4533              logger.register(es, append).add()  # initial values, not fitness values 
4534   
4535              # if es.countiter == 0 and es.opts['verb_log'] > 0 and not es.opts['verb_append']: 
4536              #    logger = CMADataLogger(es.opts['verb_filenameprefix']).register(es) 
4537              #    logger.add() 
4538              # es.writeOutput()  # initial values for sigma etc 
4539   
4540              if 1 < 3: 
4541                  if noise_handler: 
4542                      noisehandler = noise_handler 
4543                      noise_handling = True 
4544                  else: 
4545                      noisehandler = NoiseHandler(es.N, 0) 
4546                      noise_handling = False 
4547                  es.noise_handler = noisehandler 
4548   
4549              # the problem: this assumes that good solutions cannot take longer than bad ones: 
4550              # with EvalInParallel(objective_function, 2, is_feasible=opts['is_feasible']) as eval_in_parallel: 
4551              if 1 < 3: 
4552                  while not es.stop():  # iteration loop 
4553                      # X, fit = eval_in_parallel(lambda: es.ask(1)[0], es.popsize, args, repetitions=noisehandler.evaluations-1) 
4554                      X, fit = es.ask_and_eval(objective_function, args, 
4555                                               evaluations=noisehandler.evaluations, 
4556                                               aggregation=np.median)  # treats NaN with resampling 
4557                      # TODO: check args and in case use args=(noisehandler.evaluations, ) 
4558   
4559                      es.tell(X, fit)  # prepare for next iteration 
4560                      if noise_handling:  # it would be better to also use these f-evaluations in tell 
4561                          es.sigma *= noisehandler(X, fit, objective_function, es.ask, 
4562                                                   args=args)**fmin_opts['noise_change_sigma_exponent'] 
4563                          es.countevals += noisehandler.evaluations_just_done  # TODO: this is a hack, not important though 
4564                          if 1 < 3: 
4565                              es.sp.cmean *= exp(-noise_kappa_exponent * np.tanh(noisehandler.noiseS)) 
4566                              if es.sp.cmean > 1: 
4567                                  es.sp.cmean = 1 
4568   
4569                      es.disp() 
4570                      logger.add(more_data=[noisehandler.evaluations, 10**noisehandler.noiseS] if noise_handling else [], 
4571                                 modulo=1 if es.stop() and logger.modulo else None) 
4572                      if (opts['verb_log'] and opts['verb_plot'] and 
4573                            (es.countiter % max(opts['verb_plot'], opts['verb_log']) == 0 or es.stop())): 
4574                          logger.plot(324, fontsize=10) 
4575   
4576              # end while not es.stop 
4577              mean_pheno = es.gp.pheno(es.mean, into_bounds=es.boundary_handler.repair, archive=es.sent_solutions) 
4578              fmean = objective_function(mean_pheno, *args) 
4579              es.countevals += 1 
4580   
4581              es.best.update([mean_pheno], es.sent_solutions, [fmean], es.countevals) 
4582              best.update(es.best, es.sent_solutions)  # in restarted case 
4583   
4584              # final message 
4585              if opts['verb_disp']: 
4586                  es.result_pretty(irun, time.asctime(time.localtime())) 
4587   
4588              irun += 1 
4589              if irun > fmin_opts['restarts'] or 'ftarget' in es.stopdict or 'maxfevals' in es.stopdict: 
4590                  break 
4591              opts['verb_append'] = es.countevals 
4592              opts['popsize'] = fmin_opts['incpopsize'] * es.sp.popsize  # TODO: use rather options? 
4593              opts['seed'] += 1 
4594   
4595          # while irun 
4596   
4597          es.out['best'] = best  # TODO: this is a rather suboptimal type for inspection in the shell 
4598          if 1 < 3: 
4599              return es.result() + (es.stop(), es, logger) 
4600   
4601          else:  # previously: to be removed 
4602              return (best.x.copy(), best.f, es.countevals, 
4603                      dict((('stopdict', CMAStopDict(es.stopdict)) 
4604                            , ('mean', es.gp.pheno(es.mean)) 
4605                            , ('std', es.sigma * sqrt(es.dC) * es.gp.scales) 
4606                            , ('out', es.out) 
4607                            , ('opts', es.opts)  # last state of options 
4608                            , ('cma', es) 
4609                            , ('inputargs', es.inputargs) 
4610                            )) 
4611                     ) 
4612          # TODO refine output, can #args be flexible? 
4613          # is this well usable as it is now? 
4614      else:  # except KeyboardInterrupt:  # Exception, e: 
4615          if eval(str(options['verb_disp'])) > 0: 
4616              print(' in/outcomment ``raise`` in last line of cma.fmin to prevent/restore KeyboardInterrupt exception') 
4617          raise  # cave: swallowing this exception can silently mess up experiments, if ctrl-C is hit 
4618   
4619  # _____________________________________________________________________ 
4620  # _____________________________________________________________________ 
4621  # 
4622 -class BaseDataLogger(object): 
4623      """"abstract" base class for a data logger that can be used with an `OOOptimizer` 
4624   
4625      Details: attribute `modulo` is used in ``OOOptimizer.optimize`` 
4626   
4627      """ 
4628 -    def add(self, optim=None, more_data=[]): 
4629          """abstract method, add a "data point" from the state of `optim` into the 
4630          logger, the argument `optim` can be omitted if it was `register()`-ed before, 
4631          acts like an event handler""" 
4632          raise NotImplementedError() 
4633 -    def register(self, optim): 
4634          """abstract method, register an optimizer `optim`, only needed if `add()` is 
4635          called without a value for the `optim` argument""" 
4636          self.optim = optim 
4637 -    def disp(self): 
4638          """display some data trace (not implemented)""" 
4639          print('method BaseDataLogger.disp() not implemented, to be done in subclass ' + str(type(self))) 
4640 -    def plot(self): 
4641          """plot data (not implemented)""" 
4642          print('method BaseDataLogger.plot() is not implemented, to be done in subclass ' + str(type(self))) 
4643 -    def data(self): 
4644          """return logged data in a dictionary (not implemented)""" 
4645          print('method BaseDataLogger.data() is not implemented, to be done in subclass ' + str(type(self))) 
4646   
4647  # _____________________________________________________________________ 
4648  # _____________________________________________________________________ 
4649  # 
4650 -class CMADataLogger(BaseDataLogger): 
4651      """data logger for class `CMAEvolutionStrategy`. The logger is 
4652      identified by its name prefix and (over-)writes or reads according 
4653      data files. Therefore, the logger must be considered as *global* variable 
4654      with unpredictable side effects, if two loggers with the same name 
4655      and on the same working folder are used at the same time. 
4656   
4657      Examples 
4658      ======== 
4659      :: 
4660   
4661          import cma 
4662          es = cma.CMAEvolutionStrategy(...) 
4663          logger = cma.CMADataLogger().register(es) 
4664          while not es.stop(): 
4665              ... 
4666              logger.add()  # add can also take an argument 
4667   
4668          logger.plot() # or a short cut can be used: 
4669          cma.plot()  # plot data from logger with default name 
4670   
4671   
4672          logger2 = cma.CMADataLogger('just_another_filename_prefix').load() 
4673          logger2.plot() 
4674          logger2.disp() 
4675   
4676      :: 
4677   
4678          import cma 
4679          from matplotlib.pylab import * 
4680          res = cma.fmin(cma.Fcts.sphere, rand(10), 1e-0) 
4681          logger = res[-1]  # the CMADataLogger 
4682          logger.load()  # by "default" data are on disk 
4683          semilogy(logger.f[:,0], logger.f[:,5])  # plot f versus iteration, see file header 
4684          show() 
4685   
4686      Details 
4687      ======= 
4688      After loading data, the logger has the attributes `xmean`, `xrecent`, `std`, `f`, and `D`, 
4689      corresponding to ``xmean``, ``xrecentbest``, ``stddev``, ``fit``, and ``axlen`` filename 
4690      trails. 
4691   
4692      :See: `disp()`, `plot()` 
4693   
4694      """ 
4695      default_prefix = 'outcmaes' 
4696      # names = ('axlen','fit','stddev','xmean','xrecentbest') 
4697      # key_names_with_annotation = ('std', 'xmean', 'xrecent') 
4698   
4699 -    def __init__(self, name_prefix=default_prefix, modulo=1, append=False): 
4700          """initialize logging of data from a `CMAEvolutionStrategy` instance, 
4701          default ``modulo=1`` means logging with each call 
4702   
4703          """ 
4704          # super(CMAData, self).__init__({'iter':[], 'stds':[], 'D':[], 'sig':[], 'fit':[], 'xm':[]}) 
4705          # class properties: 
4706          self.file_names = ('axlen', 'fit', 'stddev', 'xmean', 'xrecentbest')  # used in load, however hard-coded in add 
4707          self.key_names = ('D', 'f', 'std', 'xmean', 'xrecent')  # used in load, however hard-coded in plot 
4708          self.key_names_with_annotation = ('std', 'xmean', 'xrecent')  # used in load 
4709          self.modulo = modulo  # allows calling with None 
4710          self.append = append 
4711          self.counter = 0  # number of calls of add, should initial value depend on `append`? 
4712          self.last_iteration = 0 
4713          self.name_prefix = name_prefix if name_prefix else CMADataLogger.default_prefix 
4714          if isinstance(self.name_prefix, CMAEvolutionStrategy): 
4715              self.name_prefix = self.name_prefix.opts.eval('verb_filenameprefix') 
4716          self.registered = False 
4717   
4718 -    def register(self, es, append=None, modulo=None): 
4719          """register a `CMAEvolutionStrategy` instance for logging, 
4720          ``append=True`` appends to previous data logged under the same name, 
4721          by default previous data are overwritten. 
4722   
4723          """ 
4724          if not isinstance(es, CMAEvolutionStrategy): 
4725              raise TypeError("only class CMAEvolutionStrategy can be registered for logging") 
4726          self.es = es 
4727          if append is not None: 
4728              self.append = append 
4729          if modulo is not None: 
4730              self.modulo = modulo 
4731          self.registered = True 
4732          return self 
4733   
4734 -    def initialize(self, modulo=None): 
4735          """reset logger, overwrite original files, `modulo`: log only every modulo call""" 
4736          if modulo is not None: 
4737              self.modulo = modulo 
4738          try: 
4739              es = self.es  # must have been registered 
4740          except AttributeError: 
4741              pass  # TODO: revise usage of es... that this can pass 
4742              raise _Error('call register() before initialize()') 
4743   
4744          self.counter = 0  # number of calls of add 
4745          self.last_iteration = 0  # some lines are only written if iteration>last_iteration 
4746   
4747          # write headers for output 
4748          fn = self.name_prefix + 'fit.dat' 
4749          strseedtime = 'seed=%d, %s' % (es.opts['seed'], time.asctime()) 
4750   
4751          try: 
4752              with open(fn, 'w') as f: 
4753                  f.write('% # columns="iteration, evaluation, sigma, axis ratio, ' + 
4754                          'bestever, best, median, worst objective function value, ' + 
4755                          'further objective values of best", ' + 
4756                          strseedtime + 
4757                          # strftime("%Y/%m/%d %H:%M:%S", localtime()) + # just asctime() would do 
4758                          '\n') 
4759          except (IOError, OSError): 
4760              print('could not open file ' + fn) 
4761   
4762          fn = self.name_prefix + 'axlen.dat' 
4763          try: 
4764              f = open(fn, 'w') 
4765              f.write('%  columns="iteration, evaluation, sigma, max axis length, ' + 
4766                      ' min axis length, all principle axes lengths ' + 
4767                      ' (sorted square roots of eigenvalues of C)", ' + 
4768                      strseedtime + 
4769                      '\n') 
4770              f.close() 
4771          except (IOError, OSError): 
4772              print('could not open file ' + fn) 
4773          finally: 
4774              f.close() 
4775          fn = self.name_prefix + 'stddev.dat' 
4776          try: 
4777              f = open(fn, 'w') 
4778              f.write('% # columns=["iteration, evaluation, sigma, void, void, ' + 
4779                      ' stds==sigma*sqrt(diag(C))", ' + 
4780                      strseedtime + 
4781                      '\n') 
4782              f.close() 
4783          except (IOError, OSError): 
4784              print('could not open file ' + fn) 
4785          finally: 
4786              f.close() 
4787   
4788          fn = self.name_prefix + 'xmean.dat' 
4789          try: 
4790              with open(fn, 'w') as f: 
4791                  f.write('% # columns="iteration, evaluation, void, void, void, xmean", ' + 
4792                          strseedtime) 
4793                  f.write(' # scaling_of_variables: ') 
4794                  if np.size(es.gp.scales) > 1: 
4795                      f.write(' '.join(map(str, es.gp.scales))) 
4796                  else: 
4797                      f.write(str(es.gp.scales)) 
4798                  f.write(', typical_x: ') 
4799                  if np.size(es.gp.typical_x) > 1: 
4800                      f.write(' '.join(map(str, es.gp.typical_x))) 
4801                  else: 
4802                      f.write(str(es.gp.typical_x)) 
4803                  f.write('\n') 
4804                  f.close() 
4805          except (IOError, OSError): 
4806              print('could not open/write file ' + fn) 
4807   
4808          fn = self.name_prefix + 'xrecentbest.dat' 
4809          try: 
4810              with open(fn, 'w') as f: 
4811                  f.write('% # iter+eval+sigma+0+fitness+xbest, ' + 
4812                          strseedtime + 
4813                          '\n') 
4814          except (IOError, OSError): 
4815              print('could not open/write file ' + fn) 
4816   
4817          return self 
4818      # end def __init__ 
4819   
4820 -    def load(self, filenameprefix=None): 
4821          """loads data from files written and return a data dictionary, *not* 
4822          a prerequisite for using `plot()` or `disp()`. 
4823   
4824          Argument `filenameprefix` is the filename prefix of data to be loaded (five files), 
4825          by default ``'outcmaes'``. 
4826   
4827          Return data dictionary with keys `xrecent`, `xmean`, `f`, `D`, `std` 
4828   
4829          """ 
4830          if not filenameprefix: 
4831              filenameprefix = self.name_prefix 
4832          for i in rglen((self.file_names)): 
4833              fn = filenameprefix + self.file_names[i] + '.dat' 
4834              try: 
4835                  self.__dict__[self.key_names[i]] = _fileToMatrix(fn) 
4836              except: 
4837                  print('WARNING: reading from file "' + fn + '" failed') 
4838              if self.key_names[i] in self.key_names_with_annotation: 
4839                  self.__dict__[self.key_names[i]].append(self.__dict__[self.key_names[i]][-1])  # copy last row to later fill in annotation position for display 
4840              self.__dict__[self.key_names[i]] = array(self.__dict__[self.key_names[i]], copy=False) 
4841          return self 
4842   
4843 -    def add(self, es=None, more_data=[], modulo=None):  # TODO: find a different way to communicate current x and f 
4844          """append some logging data from `CMAEvolutionStrategy` class instance `es`, 
4845          if ``number_of_times_called % modulo`` equals to zero, never if ``modulo==0``. 
4846   
4847          The sequence ``more_data`` must always have the same length. 
4848   
4849          When used for a different optimizer class, this function can be 
4850          (easily?) adapted by changing the assignments under INTERFACE 
4851          in the implemention. 
4852   
4853          """ 
4854          mod = modulo if modulo is not None else self.modulo 
4855          self.counter += 1 
4856          if mod == 0 or (self.counter > 3 and (self.counter - 1) % mod): 
4857              return 
4858          if es is None: 
4859              try: 
4860                  es = self.es  # must have been registered 
4861              except AttributeError : 
4862                  raise _Error('call `add` with argument `es` or ``register(es)`` before ``add()``') 
4863          elif not self.registered: 
4864              self.register(es) 
4865   
4866          if 1 < 3: 
4867              if self.counter == 1 and not self.append and self.modulo != 0: 
4868                  self.initialize()  # write file headers 
4869                  self.counter = 1 
4870   
4871          # --- INTERFACE, can be changed if necessary --- 
4872          if not isinstance(es, CMAEvolutionStrategy):  # not necessary 
4873              print('WARNING: <type \'CMAEvolutionStrategy\'> expected, found ' 
4874                              + str(type(es)) + ' in method CMADataLogger.add') 
4875          evals = es.countevals 
4876          iteration = es.countiter 
4877          sigma = es.sigma 
4878          axratio = es.D.max() / es.D.min() 
4879          xmean = es.mean  # TODO: should be optionally phenotype? 
4880          fmean_noise_free = es.fmean_noise_free 
4881          fmean = es.fmean 
4882          try: 
4883              besteverf = es.best.f 
4884              bestf = es.fit.fit[0] 
4885              worstf = es.fit.fit[-1] 
4886              medianf = es.fit.fit[es.sp.popsize // 2] 
4887          except: 
4888              if iteration > 0:  # first call without f-values is OK 
4889                  raise 
4890          try: 
4891              xrecent = es.best.last.x 
4892          except: 
4893              xrecent = None 
4894          maxD = es.D.max() 
4895          minD = es.D.min() 
4896          diagD = es.D 
4897          diagC = es.sigma * es.sigma_vec * sqrt(es.dC) 
4898          more_to_write = es.more_to_write 
4899          es.more_to_write = [] 
4900          # --- end interface --- 
4901   
4902          try: 
4903              # fit 
4904              if iteration > self.last_iteration: 
4905                  fn = self.name_prefix + 'fit.dat' 
4906                  with open(fn, 'a') as f: 
4907                      f.write(str(iteration) + ' ' 
4908                              + str(evals) + ' ' 
4909                              + str(sigma) + ' ' 
4910                              + str(axratio) + ' ' 
4911                              + str(besteverf) + ' ' 
4912                              + '%.16e' % bestf + ' ' 
4913                              + str(medianf) + ' ' 
4914                              + str(worstf) + ' ' 
4915                              # + str(es.sp.popsize) + ' ' 
4916                              # + str(10**es.noiseS) + ' ' 
4917                              # + str(es.sp.cmean) + ' ' 
4918                              + ' '.join(str(i) for i in more_to_write) + ' ' 
4919                              + ' '.join(str(i) for i in more_data) + ' ' 
4920                              + '\n') 
4921              # axlen 
4922              fn = self.name_prefix + 'axlen.dat' 
4923              with open(fn, 'a') as f:  # does not rely on reference counting 
4924                  f.write(str(iteration) + ' ' 
4925                          + str(evals) + ' ' 
4926                          + str(sigma) + ' ' 
4927                          + str(maxD) + ' ' 
4928                          + str(minD) + ' ' 
4929                          + ' '.join(map(str, diagD)) 
4930                          + '\n') 
4931              # stddev 
4932              fn = self.name_prefix + 'stddev.dat' 
4933              with open(fn, 'a') as f: 
4934                  f.write(str(iteration) + ' ' 
4935                          + str(evals) + ' ' 
4936                          + str(sigma) + ' ' 
4937                          + '0 0 ' 
4938                          + ' '.join(map(str, diagC)) 
4939                          + '\n') 
4940              # xmean 
4941              fn = self.name_prefix + 'xmean.dat' 
4942              with open(fn, 'a') as f: 
4943                  f.write(str(iteration) + ' ' 
4944                          + str(evals) + ' ' 
4945                          # + str(sigma) + ' ' 
4946                          + '0 ' 
4947                          + str(fmean_noise_free) + ' ' 
4948                          + str(fmean) + ' '  # TODO: this does not make sense 
4949                          # TODO should be optional the phenotyp? 
4950                          + ' '.join(map(str, xmean)) 
4951                          + '\n') 
4952              # xrecent 
4953              fn = self.name_prefix + 'xrecentbest.dat' 
4954              if iteration > 0 and xrecent is not None: 
4955                  with open(fn, 'a') as f: 
4956                      f.write(str(iteration) + ' ' 
4957                              + str(evals) + ' ' 
4958                              + str(sigma) + ' ' 
4959                              + '0 ' 
4960                              + str(bestf) + ' ' 
4961                              + ' '.join(map(str, xrecent)) 
4962                              + '\n') 
4963   
4964          except (IOError, OSError): 
4965              if iteration <= 1: 
4966                  print('could not open/write file') 
4967          self.last_iteration = iteration 
4968   
4969 -    def closefig(self): 
4970          pyplot.close(self.fighandle) 
4971   
4972 -    def save(self, nameprefix, switch=False): 
4973          """saves logger data to a different set of files, for 
4974          ``switch=True`` also the loggers name prefix is switched to 
4975          the new value 
4976   
4977          """ 
4978          if not nameprefix or not isinstance(nameprefix, basestring): 
4979              raise _Error('filename prefix must be a nonempty string') 
4980   
4981          if nameprefix == self.default_prefix: 
4982              raise _Error('cannot save to default name "' + nameprefix + '...", chose another name') 
4983   
4984          if nameprefix == self.name_prefix: 
4985              return 
4986   
4987          for name in CMADataLogger.names: 
4988              open(nameprefix + name + '.dat', 'w').write(open(self.name_prefix + name + '.dat').read()) 
4989   
4990          if switch: 
4991              self.name_prefix = nameprefix 
4992   
4993 -    def plot(self, fig=None, iabscissa=1, iteridx=None, 
4994               plot_mean=False,  # was: plot_mean=True 
4995               foffset=1e-19, x_opt=None, fontsize=10): 
4996          """ 
4997          plot data from a `CMADataLogger` (using the files written by the logger). 
4998   
4999          Arguments 
5000          --------- 
5001              `fig` 
5002                  figure number, by default 325 
5003              `iabscissa` 
5004                  ``0==plot`` versus iteration count, 
5005                  ``1==plot`` versus function evaluation number 
5006              `iteridx` 
5007                  iteration indices to plot 
5008   
5009          Return `CMADataLogger` itself. 
5010   
5011          Examples 
5012          -------- 
5013          :: 
5014   
5015              import cma 
5016              logger = cma.CMADataLogger()  # with default name 
5017              # try to plot the "default logging" data (e.g. 
5018              #   from previous fmin calls, which is essentially what 
5019              #   also cma.plot() does) 
5020              logger.plot() 
5021              cma.savefig('fig325.png')  # save current figure 
5022              logger.closefig() 
5023   
5024          Dependencies: matlabplotlib/pyplot. 
5025   
5026          """ 
5027          dat = self.load(self.name_prefix) 
5028          try: 
5029              # pyplot: prodedural interface for matplotlib 
5030              from  matplotlib.pyplot import figure, subplot, semilogy, hold, plot, grid, \ 
5031                   axis, title, text, xlabel, isinteractive, gcf 
5032   
5033          except ImportError: 
5034              ImportError('could not find matplotlib.pyplot module, function plot() is not available') 
5035              return 
5036   
5037          if fontsize and pyplot.rcParams['font.size'] != fontsize: 
5038              print('global variable pyplot.rcParams[\'font.size\'] set (from ' + 
5039                    str(pyplot.rcParams['font.size']) + ') to ' + str(fontsize)) 
5040              pyplot.rcParams['font.size'] = fontsize  # subtracted in the end, but return can happen inbetween 
5041   
5042          if fig: 
5043              figure(fig) 
5044          else: 
5045              figure(325) 
5046              # show()  # should not be necessary 
5047          self.fighandle = gcf()  # fighandle.number 
5048   
5049          if iabscissa not in (0, 1): 
5050              iabscissa = 1 
5051   
5052          # interactive_status = matplotlib.is_interactive() 
5053          pyplot.ioff()  # prevents immediate drawing, much faster 
5054   
5055          dat.x = dat.xmean  # this is the genotyp 
5056          if not plot_mean: 
5057              if len(dat.x) < 2: 
5058                  print('not enough data to plot recent x') 
5059              else: 
5060                  dat.x = dat.xrecent 
5061   
5062          if iteridx is not None: 
5063              dat.f = dat.f[np.where([x in iteridx for x in dat.f[:, 0]])[0], :] 
5064              dat.D = dat.D[np.where([x in iteridx for x in dat.D[:, 0]])[0], :] 
5065              iteridx.append(dat.x[-1, 1])  # last entry is artificial 
5066              dat.x = dat.x[np.where([x in iteridx for x in dat.x[:, 0]])[0], :] 
5067              dat.std = dat.std[np.where([x in iteridx for x in dat.std[:, 0]])[0], :] 
5068   
5069          if iabscissa == 0: 
5070              xlab = 'iterations' 
5071          elif iabscissa == 1: 
5072              xlab = 'function evaluations' 
5073   
5074          # use fake last entry in x and std for line extension-annotation 
5075          if dat.x.shape[1] < 100: 
5076              minxend = int(1.06 * dat.x[-2, iabscissa]) 
5077              # write y-values for individual annotation into dat.x 
5078              dat.x[-1, iabscissa] = minxend  # TODO: should be ax[1] 
5079              idx = np.argsort(dat.x[-2, 5:]) 
5080              idx2 = np.argsort(idx) 
5081              if x_opt is None: 
5082                  dat.x[-1, 5 + idx] = np.linspace(np.min(dat.x[:, 5:]), 
5083                              np.max(dat.x[:, 5:]), dat.x.shape[1] - 5) 
5084              else: 
5085                  dat.x[-1, 5 + idx] = np.logspace(np.log10(np.min(abs(dat.x[:, 5:]))), 
5086                              np.log10(np.max(abs(dat.x[:, 5:]))), dat.x.shape[1] - 5) 
5087          else: 
5088              minxend = 0 
5089   
5090          if len(dat.f) == 0: 
5091              print('nothing to plot') 
5092              return 
5093   
5094          # not in use anymore, see formatter above 
5095          # xticklocs = np.arange(5) * np.round(minxend/4., -int(np.log10(minxend/4.))) 
5096   
5097          # dfit(dfit<1e-98) = NaN; 
5098   
5099          # TODO: if abscissa==0 plot in chunks, ie loop over subsets where dat.f[:,0]==countiter is monotonous 
5100   
5101          subplot(2, 2, 1) 
5102          self.plotdivers(dat, iabscissa, foffset) 
5103   
5104          # TODO: modularize also the remaining subplots 
5105          subplot(2, 2, 2) 
5106          hold(False) 
5107          if x_opt is not None:  # TODO: differentate neg and pos? 
5108              semilogy(dat.x[:, iabscissa], abs(dat.x[:, 5:]) - x_opt, '-') 
5109          else: 
5110              plot(dat.x[:, iabscissa], dat.x[:, 5:], '-') 
5111          hold(True) 
5112          grid(True) 
5113          ax = array(axis()) 
5114          # ax[1] = max(minxend, ax[1]) 
5115          axis(ax) 
5116          ax[1] -= 1e-6 
5117          if dat.x.shape[1] < 100: 
5118              yy = np.linspace(ax[2] + 1e-6, ax[3] - 1e-6, dat.x.shape[1] - 5) 
5119              # yyl = np.sort(dat.x[-1,5:]) 
5120              idx = np.argsort(dat.x[-1, 5:]) 
5121              idx2 = np.argsort(idx) 
5122              if x_opt is not None: 
5123                  semilogy([dat.x[-1, iabscissa], ax[1]], [abs(dat.x[-1, 5:]), yy[idx2]], 'k-')  # line from last data point 
5124                  semilogy(np.dot(dat.x[-2, iabscissa], [1, 1]), array([ax[2] + 1e-6, ax[3] - 1e-6]), 'k-') 
5125              else: 
5126                  # plot([dat.x[-1, iabscissa], ax[1]], [dat.x[-1,5:], yy[idx2]], 'k-') # line from last data point 
5127                  plot(np.dot(dat.x[-2, iabscissa], [1, 1]), array([ax[2] + 1e-6, ax[3] - 1e-6]), 'k-') 
5128              # plot(array([dat.x[-1, iabscissa], ax[1]]), 
5129              #      reshape(array([dat.x[-1,5:], yy[idx2]]).flatten(), (2,4)), '-k') 
5130              for i in range(len(idx)): 
5131                  # TODOqqq: annotate phenotypic value!? 
5132                  # text(ax[1], yy[i], 'x(' + str(idx[i]) + ')=' + str(dat.x[-2,5+idx[i]])) 
5133                  text(dat.x[-1, iabscissa], dat.x[-1, 5 + i], 'x(' + str(i) + ')=' + str(dat.x[-2, 5 + i])) 
5134   
5135          i = 2  # find smallest i where iteration count differs (in case the same row appears twice) 
5136          while i < len(dat.f) and dat.f[-i][0] == dat.f[-1][0]: 
5137              i += 1 
5138          title('Object Variables (' + ('mean' if plot_mean else 'curr best') + 
5139                  ', ' + str(dat.x.shape[1] - 5) + '-D, popsize~' + 
5140                  (str(int((dat.f[-1][1] - dat.f[-i][1]) / (dat.f[-1][0] - dat.f[-i][0]))) 
5141                      if len(dat.f.T[0]) > 1 and dat.f[-1][0] > dat.f[-i][0] else 'NA') 
5142                  + ')') 
5143          # pyplot.xticks(xticklocs) 
5144   
5145          # Scaling 
5146          subplot(2, 2, 3) 
5147          hold(False) 
5148          semilogy(dat.D[:, iabscissa], dat.D[:, 5:], '-b') 
5149          hold(True) 
5150          grid(True) 
5151          ax = array(axis()) 
5152          # ax[1] = max(minxend, ax[1]) 
5153          axis(ax) 
5154          title('Scaling (All Main Axes)') 
5155          # pyplot.xticks(xticklocs) 
5156          xlabel(xlab) 
5157   
5158          # standard deviations 
5159          subplot(2, 2, 4) 
5160          hold(False) 
5161          # remove sigma from stds (graphs become much better readible) 
5162          dat.std[:, 5:] = np.transpose(dat.std[:, 5:].T / dat.std[:, 2].T) 
5163          # ax = array(axis()) 
5164          # ax[1] = max(minxend, ax[1]) 
5165          # axis(ax) 
5166          if 1 < 2 and dat.std.shape[1] < 100: 
5167              # use fake last entry in x and std for line extension-annotation 
5168              minxend = int(1.06 * dat.x[-2, iabscissa]) 
5169              dat.std[-1, iabscissa] = minxend  # TODO: should be ax[1] 
5170              idx = np.argsort(dat.std[-2, 5:]) 
5171              idx2 = np.argsort(idx) 
5172              dat.std[-1, 5 + idx] = np.logspace(np.log10(np.min(dat.std[:, 5:])), 
5173                              np.log10(np.max(dat.std[:, 5:])), dat.std.shape[1] - 5) 
5174   
5175              dat.std[-1, iabscissa] = minxend  # TODO: should be ax[1] 
5176              yy = np.logspace(np.log10(ax[2]), np.log10(ax[3]), dat.std.shape[1] - 5) 
5177              # yyl = np.sort(dat.std[-1,5:]) 
5178              idx = np.argsort(dat.std[-1, 5:]) 
5179              idx2 = np.argsort(idx) 
5180              # plot(np.dot(dat.std[-2, iabscissa],[1,1]), array([ax[2]+1e-6, ax[3]-1e-6]), 'k-') # vertical separator 
5181              # vertical separator 
5182              plot(np.dot(dat.std[-2, iabscissa], [1, 1]), array([np.min(dat.std[-2, 5:]), np.max(dat.std[-2, 5:])]), 'k-') 
5183              hold(True) 
5184              # plot([dat.std[-1, iabscissa], ax[1]], [dat.std[-1,5:], yy[idx2]], 'k-') # line from last data point 
5185              for i in rglen((idx)): 
5186                  # text(ax[1], yy[i], ' '+str(idx[i])) 
5187                  text(dat.std[-1, iabscissa], dat.std[-1, 5 + i], ' ' + str(i)) 
5188          semilogy(dat.std[:, iabscissa], dat.std[:, 5:], '-') 
5189          grid(True) 
5190          title('Standard Deviations in All Coordinates') 
5191          # pyplot.xticks(xticklocs) 
5192          xlabel(xlab) 
5193          pyplot.ion() 
5194          pyplot.draw()  # update "screen" 
5195          pyplot.show()  # show figure 
5196          # matplotlib.interactive(interactive_status) 
5197   
5198          return self 
5199   
5200      # ____________________________________________________________ 
5201      # ____________________________________________________________ 
5202      # 
5203      @staticmethod 
5204 -    def plotdivers(dat, iabscissa, foffset): 
5205          """helper function for `plot()` that plots all what is 
5206          in the upper left subplot like fitness, sigma, etc. 
5207   
5208          Arguments 
5209          --------- 
5210              `iabscissa` in ``(0,1)`` 
5211                  0==versus fevals, 1==versus iteration 
5212              `foffset` 
5213                  offset to fitness for log-plot 
5214   
5215           :See: `plot()` 
5216   
5217          """ 
5218          from matplotlib.pyplot import semilogy, hold, grid, \ 
5219              axis, title, text 
5220          fontsize = pyplot.rcParams['font.size'] 
5221   
5222          # interactive_status = matplotlib.is_interactive() 
5223          pyplot.ioff()  # prevents immediate drawing 
5224          hold(False) 
5225   
5226          dfit = dat.f[:, 5] - min(dat.f[:, 5]) 
5227          dfit[dfit < 1e-98] = np.NaN 
5228   
5229          if dat.f.shape[1] > 7: 
5230              # semilogy(dat.f[:, iabscissa], abs(dat.f[:,[6, 7, 10, 12]])+foffset,'-k') 
5231              semilogy(dat.f[:, iabscissa], abs(dat.f[:, [6, 7]]) + foffset, '-k') 
5232              hold(True) 
5233   
5234          # (larger indices): additional fitness data, for example constraints values 
5235          if dat.f.shape[1] > 8: 
5236              # dd = abs(dat.f[:,7:]) + 10*foffset 
5237              # dd = np.where(dat.f[:,7:]==0, np.NaN, dd) # cannot be 
5238              semilogy(dat.f[:, iabscissa], np.abs(dat.f[:, 8:]) + 10 * foffset, 'm') 
5239              hold(True) 
5240   
5241          idx = np.where(dat.f[:, 5] > 1e-98)[0]  # positive values 
5242          semilogy(dat.f[idx, iabscissa], dat.f[idx, 5] + foffset, '.b') 
5243          hold(True) 
5244          grid(True) 
5245   
5246          idx = np.where(dat.f[:, 5] < -1e-98)[0]  # negative values 
5247          semilogy(dat.f[idx, iabscissa], abs(dat.f[idx, 5]) + foffset, '.r') 
5248   
5249          semilogy(dat.f[:, iabscissa], abs(dat.f[:, 5]) + foffset, '-b') 
5250          idx = np.isfinite(dfit) 
5251          semilogy(dat.f[:, iabscissa][idx], dfit[idx], '-c') 
5252   
5253          # overall minimum 
5254          i = np.argmin(dat.f[:, 5]) 
5255          semilogy(dat.f[i, iabscissa] * np.ones(2), dat.f[i, 5] * np.ones(2), 'rd') 
5256          # semilogy(dat.f[-1, iabscissa]*np.ones(2), dat.f[-1,4]*np.ones(2), 'rd') 
5257   
5258          # AR and sigma 
5259          semilogy(dat.f[:, iabscissa], dat.f[:, 3], '-r')  # AR 
5260          semilogy(dat.f[:, iabscissa], dat.f[:, 2], '-g')  # sigma 
5261          semilogy(dat.std[:-1, iabscissa], np.vstack([list(map(max, dat.std[:-1, 5:])), list(map(min, dat.std[:-1, 5:]))]).T, 
5262                       '-m', linewidth=2) 
5263          text(dat.std[-2, iabscissa], max(dat.std[-2, 5:]), 'max std', fontsize=fontsize) 
5264          text(dat.std[-2, iabscissa], min(dat.std[-2, 5:]), 'min std', fontsize=fontsize) 
5265          ax = array(axis()) 
5266          # ax[1] = max(minxend, ax[1]) 
5267          axis(ax) 
5268          text(ax[0] + 0.01, ax[2],  # 10**(log10(ax[2])+0.05*(log10(ax[3])-log10(ax[2]))), 
5269               '.f_recent=' + repr(dat.f[-1, 5])) 
5270   
5271          # title('abs(f) (blue), f-min(f) (cyan), Sigma (green), Axis Ratio (red)') 
5272          title('blue:abs(f), cyan:f-min(f), green:sigma, red:axis ratio', fontsize=fontsize - 1) 
5273          # pyplot.xticks(xticklocs) 
5274          pyplot.ion() 
5275          pyplot.draw()  # update "screen" 
5276          pyplot.show()  # show figure 
5277          # matplotlib.interactive(interactive_status) 
5278   
5279   
5280 -    def downsampling(self, factor=10, first=3, switch=True, verbose=True): 
5281          """ 
5282          rude downsampling of a `CMADataLogger` data file by `factor`, keeping 
5283          also the first `first` entries. This function is a stump and subject 
5284          to future changes. Return self. 
5285   
5286          Arguments 
5287          --------- 
5288             - `factor` -- downsampling factor 
5289             - `first` -- keep first `first` entries 
5290             - `switch` -- switch the new logger to the downsampled logger original_name+'down' 
5291   
5292          Details 
5293          ------- 
5294          ``self.name_prefix+'down'`` files are written 
5295   
5296          Example 
5297          ------- 
5298          :: 
5299   
5300              import cma 
5301              cma.downsampling()  # takes outcmaes* files 
5302              cma.plot('outcmaesdown') 
5303   
5304          """ 
5305          newprefix = self.name_prefix + 'down' 
5306          for name in self.file_names: 
5307              f = open(newprefix + name + '.dat', 'w') 
5308              iline = 0 
5309              cwritten = 0 
5310              for line in open(self.name_prefix + name + '.dat'): 
5311                  if iline < first or iline % factor == 0: 
5312                      f.write(line) 
5313                      cwritten += 1 
5314                  iline += 1 
5315              f.close() 
5316              if verbose and iline > first: 
5317                  print('%d' % (cwritten) + ' lines written in ' + newprefix + name + '.dat') 
5318          if switch: 
5319              self.name_prefix += 'down' 
5320          return self 
5321   
5322      # ____________________________________________________________ 
5323      # ____________________________________________________________ 
5324      # 
5325 -    def disp(self, idx=100):  # r_[0:5,1e2:1e9:1e2,-10:0]): 
5326          """displays selected data from (files written by) the class `CMADataLogger`. 
5327   
5328          Arguments 
5329          --------- 
5330             `idx` 
5331                 indices corresponding to rows in the data file; 
5332                 if idx is a scalar (int), the first two, then every idx-th, 
5333                 and the last three rows are displayed. Too large index values are removed. 
5334   
5335          Example 
5336          ------- 
5337          >>> import cma, numpy as np 
5338          >>> res = cma.fmin(cma.fcts.elli, 7 * [0.1], 1, {'verb_disp':1e9})  # generate data 
5339          >>> assert res[1] < 1e-9 
5340          >>> assert res[2] < 4400 
5341          >>> l = cma.CMADataLogger()  # == res[-1], logger with default name, "points to" above data 
5342          >>> l.disp([0,-1])  # first and last 
5343          >>> l.disp(20)  # some first/last and every 20-th line 
5344          >>> l.disp(np.r_[0:999999:100, -1]) # every 100-th and last 
5345          >>> l.disp(np.r_[0, -10:0]) # first and ten last 
5346          >>> cma.disp(l.name_prefix, np.r_[0::100, -10:])  # the same as l.disp(...) 
5347   
5348          Details 
5349          ------- 
5350          The data line with the best f-value is displayed as last line. 
5351   
5352          :See: `disp()` 
5353   
5354          """ 
5355   
5356          filenameprefix = self.name_prefix 
5357   
5358          def printdatarow(dat, iteration): 
5359              """print data of iteration i""" 
5360              i = np.where(dat.f[:, 0] == iteration)[0][0] 
5361              j = np.where(dat.std[:, 0] == iteration)[0][0] 
5362              print('%5d' % (int(dat.f[i, 0])) + ' %6d' % (int(dat.f[i, 1])) + ' %.14e' % (dat.f[i, 5]) + 
5363                    ' %5.1e' % (dat.f[i, 3]) + 
5364                    ' %6.2e' % (max(dat.std[j, 5:])) + ' %6.2e' % min(dat.std[j, 5:])) 
5365   
5366          dat = CMADataLogger(filenameprefix).load() 
5367          ndata = dat.f.shape[0] 
5368   
5369          # map index to iteration number, is difficult if not all iteration numbers exist 
5370          # idx = idx[np.where(map(lambda x: x in dat.f[:,0], idx))[0]] # TODO: takes pretty long 
5371          # otherwise: 
5372          if idx is None: 
5373              idx = 100 
5374          if np.isscalar(idx): 
5375              # idx = np.arange(0, ndata, idx) 
5376              if idx: 
5377                  idx = np.r_[0, 1, idx:ndata - 3:idx, -3:0] 
5378              else: 
5379                  idx = np.r_[0, 1, -3:0] 
5380   
5381          idx = array(idx) 
5382          idx = idx[idx < ndata] 
5383          idx = idx[-idx <= ndata] 
5384          iters = dat.f[idx, 0] 
5385          idxbest = np.argmin(dat.f[:, 5]) 
5386          iterbest = dat.f[idxbest, 0] 
5387   
5388          if len(iters) == 1: 
5389              printdatarow(dat, iters[0]) 
5390          else: 
5391              self.disp_header() 
5392              for i in iters: 
5393                  printdatarow(dat, i) 
5394              self.disp_header() 
5395              printdatarow(dat, iterbest) 
5396          sys.stdout.flush() 
5397 -    def disp_header(self): 
5398          heading = 'Iterat Nfevals  function value    axis ratio maxstd  minstd' 
5399          print(heading) 
5400   
5401  # end class CMADataLogger 
5402   
5403  # ____________________________________________________________ 
5404  # ____________________________________________________________ 
5405  # 
5406  # _____________________________________________________________________ 
5407  # _____________________________________________________________________ 
5408  # 
5409 -class DEAPCMADataLogger(BaseDataLogger): 
5410      """data logger for class `deap.cma.Strategy` of the DEAP library. 
5411      Documentation is not entirely adapted to the deap case. 
5412   
5413      The logger is identified by its name prefix and writes or reads according 
5414      data files. 
5415   
5416      Examples 
5417      ======== 
5418      :: 
5419   
5420          import cma_logger 
5421          es = deap.cma.Strategy(...) 
5422          data = cma.DEAPCMADataLogger().register(es) 
5423          while not es.stop(): 
5424              ... 
5425              data.add(fitness_values)  # add can also take `es` as additional argument 
5426   
5427          data.plot() # or a short cut can be used: 
5428          cma.plot()  # plot data from logger with default name 
5429   
5430   
5431          data2 = cma.DEAPCMADataLogger(another_filename_prefix).load() 
5432          data2.plot() 
5433          data2.disp() 
5434   
5435      :: 
5436   
5437          import cma 
5438          from pyplot import * 
5439          res = cma.fmin(cma.Fcts.sphere, rand(10), 1e-0) 
5440          dat = res[-1]  # the CMADataLogger 
5441          dat.load()  # by "default" data are on disk 
5442          semilogy(dat.f[:,0], dat.f[:,5])  # plot f versus iteration, see file header 
5443          show() 
5444   
5445      Details 
5446      ======= 
5447      After loading data, the logger has the attributes `xmean`, `xrecent`, `std`, `f`, and `D`, 
5448      corresponding to xmean, xrecentbest, stddev, fit, and axlen filename trails. 
5449   
5450      :See: `disp()`, `plot()` 
5451   
5452      """ 
5453      default_prefix = 'outcmaes' 
5454      names = ('axlen', 'fit', 'stddev', 'xmean')  # ,'xrecentbest') 
5455      key_names_with_annotation = ('std', 'xmean') 
5456   
5457 -    def __init__(self, name_prefix=default_prefix, modulo=1, append=False): 
5458          """initialize logging of data from a `CMAEvolutionStrategy` instance, 
5459          default modulo expands to 1 == log with each call 
5460   
5461          """ 
5462          # super(CMAData, self).__init__({'iter':[], 'stds':[], 'D':[], 'sig':[], 'fit':[], 'xm':[]}) 
5463          # class properties: 
5464          self.counter = 0  # number of calls of add 
5465          self.best_fitness = np.inf 
5466          self.modulo = modulo  # allows calling with None 
5467          self.append = append 
5468          self.name_prefix = name_prefix if name_prefix else CMADataLogger.default_prefix 
5469          if isinstance(self.name_prefix, CMAEvolutionStrategy): 
5470              self.name_prefix = self.name_prefix.opts.eval('verb_filenameprefix') 
5471          self.registered = False 
5472   
5473 -    def register(self, es, append=None, modulo=None): 
5474          """register a `CMAEvolutionStrategy` instance for logging, 
5475          ``append=True`` appends to previous data logged under the same name, 
5476          by default previous data are overwritten. 
5477   
5478          """ 
5479          self.es = es 
5480          if append is not None: 
5481              self.append = append 
5482          if modulo is not None: 
5483              self.modulo = modulo 
5484          if not self.append and self.modulo != 0: 
5485              self.initialize()  # write file headers 
5486          self.registered = True 
5487          return self 
5488   
5489 -    def initialize(self, modulo=None): 
5490          """reset logger, overwrite original files, `modulo`: log only every modulo call""" 
5491          if modulo is not None: 
5492              self.modulo = modulo 
5493          try: 
5494              es = self.es  # must have been registered 
5495          except AttributeError: 
5496              pass  # TODO: revise usage of es... that this can pass 
5497              raise _Error('call register() before initialize()') 
5498   
5499          # write headers for output 
5500          fn = self.name_prefix + 'fit.dat' 
5501          if 1 < 3: 
5502              strseedtime = 'seed=unkown, %s' % (time.asctime()) 
5503   
5504          try: 
5505              with open(fn, 'w') as f: 
5506                  f.write('% # columns="iteration, evaluation, sigma, axis ratio, ' + 
5507                          'bestever, best, median, worst objective function value, ' + 
5508                          'further objective values of best", ' + 
5509                          strseedtime + 
5510                          # strftime("%Y/%m/%d %H:%M:%S", localtime()) + # just asctime() would do 
5511                          '\n') 
5512          except (IOError, OSError): 
5513              print('could not open file ' + fn) 
5514   
5515          fn = self.name_prefix + 'axlen.dat' 
5516          try: 
5517              f = open(fn, 'w') 
5518              f.write('%  columns="iteration, evaluation, sigma, max axis length, ' + 
5519                      ' min axis length, all principle axes lengths ' + 
5520                      ' (sorted square roots of eigenvalues of C)", ' + 
5521                      strseedtime + 
5522                      '\n') 
5523              f.close() 
5524          except (IOError, OSError): 
5525              print('could not open file ' + fn) 
5526          finally: 
5527              f.close() 
5528          fn = self.name_prefix + 'stddev.dat' 
5529          try: 
5530              f = open(fn, 'w') 
5531              f.write('% # columns=["iteration, evaluation, sigma, void, void, ' + 
5532                      ' stds==sigma*sqrt(diag(C))", ' + 
5533                      strseedtime + 
5534                      '\n') 
5535              f.close() 
5536          except (IOError, OSError): 
5537              print('could not open file ' + fn) 
5538          finally: 
5539              f.close() 
5540   
5541          fn = self.name_prefix + 'xmean.dat' 
5542          try: 
5543              with open(fn, 'w') as f: 
5544                  f.write('% # columns="iteration, evaluation, void, void, void, xmean", ' + 
5545                          strseedtime) 
5546                  f.write('\n') 
5547                  f.close() 
5548          except (IOError, OSError): 
5549              print('could not open/write file ' + fn) 
5550   
5551          return self 
5552      # end def __init__ 
5553   
5554 -    def load(self, filenameprefix=None): 
5555          """loads data from files written and return a data dictionary, *not* 
5556          a prerequisite for using `plot()` or `disp()`. 
5557   
5558          Argument `filenameprefix` is the filename prefix of data to be loaded (five files), 
5559          by default ``'outcmaes'``. 
5560   
5561          Return data dictionary with keys `xrecent`, `xmean`, `f`, `D`, `std` 
5562   
5563          """ 
5564          if not filenameprefix: 
5565              filenameprefix = self.name_prefix 
5566          dat = self  # historical 
5567          # dat.xrecent = _fileToMatrix(filenameprefix + 'xrecentbest.dat') 
5568          dat.xmean = _fileToMatrix(filenameprefix + 'xmean.dat') 
5569          dat.std = _fileToMatrix(filenameprefix + 'stddev' + '.dat') 
5570          # a hack to later write something into the last entry 
5571          for key in ['xmean', 'std']:  # 'xrecent', 
5572              dat.__dict__[key].append(dat.__dict__[key][-1])  # copy last row to later fill in annotation position for display 
5573              dat.__dict__[key] = array(dat.__dict__[key], copy=False) 
5574          dat.f = array(_fileToMatrix(filenameprefix + 'fit.dat')) 
5575          dat.D = array(_fileToMatrix(filenameprefix + 'axlen' + '.dat')) 
5576          return dat 
5577   
5578   
5579 -    def add(self, fitness_values, es=None, more_data=[], modulo=None):  # TODO: find a different way to communicate current x and f 
5580          """append some logging data from `CMAEvolutionStrategy` class instance `es`, 
5581          if ``number_of_times_called % modulo`` equals to zero, never if ``modulo==0``. 
5582   
5583          The sequence ``more_data`` must always have the same length. 
5584   
5585          """ 
5586          self.counter += 1 
5587          fitness_values = np.sort(fitness_values) 
5588          if fitness_values[0] < self.best_fitness: 
5589              self.best_fitness = fitness_values[0] 
5590          mod = modulo if modulo is not None else self.modulo 
5591          if mod == 0 or (self.counter > 3 and self.counter % mod): 
5592              return 
5593          if es is None: 
5594              try: 
5595                  es = self.es  # must have been registered 
5596              except AttributeError : 
5597                  raise _Error('call register() before add() or add(es)') 
5598          elif not self.registered: 
5599              self.register(es) 
5600   
5601          try: 
5602              # fit 
5603              if es.update_count > 0: 
5604                  # fit = es.fit.fit[0]  # TODO: where do we get the fitness from? 
5605                  fn = self.name_prefix + 'fit.dat' 
5606                  with open(fn, 'a') as f: 
5607                      f.write(str(es.update_count) + ' ' 
5608                              + str(es.update_count * es.lambda_) + ' ' 
5609                              + str(es.sigma) + ' ' 
5610                              + str(es.diagD[-1] / es.diagD[0]) + ' ' 
5611                              + str(self.best_fitness) + ' ' 
5612                              + '%.16e' % fitness_values[0] + ' ' 
5613                              + str(fitness_values[es.lambda_ // 2]) + ' ' 
5614                              + str(fitness_values[-1]) + ' ' 
5615                              # + str(es.sp.popsize) + ' ' 
5616                              # + str(10**es.noiseS) + ' ' 
5617                              # + str(es.sp.cmean) + ' ' 
5618                              # + ' '.join(str(i) for i in es.more_to_write) 
5619                              + ' '.join(str(i) for i in more_data) 
5620                              + '\n') 
5621                      # es.more_to_write = [] 
5622              # axlen 
5623              fn = self.name_prefix + 'axlen.dat' 
5624              with open(fn, 'a') as f:  # does not rely on reference counting 
5625                  f.write(str(es.update_count) + ' ' 
5626                          + str(es.update_count * es.lambda_) + ' ' 
5627                          + str(es.sigma) + ' ' 
5628                          + str(es.diagD[-1]) + ' ' 
5629                          + str(es.diagD[0]) + ' ' 
5630                          + ' '.join(map(str, es.diagD)) 
5631                          + '\n') 
5632              # stddev 
5633              fn = self.name_prefix + 'stddev.dat' 
5634              with open(fn, 'a') as f: 
5635                  f.write(str(es.update_count) + ' ' 
5636                          + str(es.update_count * es.lambda_) + ' ' 
5637                          + str(es.sigma) + ' ' 
5638                          + '0 0 ' 
5639                          + ' '.join(map(str, es.sigma * np.sqrt([es.C[i][i] for i in xrange(es.dim)]))) 
5640                          + '\n') 
5641              # xmean 
5642              fn = self.name_prefix + 'xmean.dat' 
5643              with open(fn, 'a') as f: 
5644                  if es.update_count < 1: 
5645                      f.write('0 0 0 0 0 ' 
5646                              + ' '.join(map(str, 
5647                                                # TODO should be optional the phenotyp? 
5648                                                # es.x0 
5649                                                es.mean)) 
5650                              + '\n') 
5651                  else: 
5652                      f.write(str(es.update_count) + ' ' 
5653                              + str(es.update_count * es.lambda_) + ' ' 
5654                              # + str(es.sigma) + ' ' 
5655                              + '0 0 0 ' 
5656                              # + str(es.fmean_noise_free) + ' ' 
5657                              # + str(es.fmean) + ' '  # TODO: this does not make sense 
5658                              # TODO should be optional the phenotyp? 
5659                              + ' '.join(map(str, es.centroid)) 
5660                              + '\n') 
5661              # xrecent 
5662          except (IOError, OSError): 
5663              if es.countiter == 1: 
5664                  print('could not open/write file') 
5665   
5666 -    def closefig(self): 
5667          pyplot.close(self.fighandle) 
5668   
5669 -    def save(self, nameprefix, switch=False): 
5670          """saves logger data to a different set of files, for 
5671          ``switch=True`` also the loggers name prefix is switched to 
5672          the new value 
5673   
5674          """ 
5675          if not nameprefix or not isinstance(nameprefix, basestring): 
5676              raise _Error('filename prefix must be a nonempty string') 
5677   
5678          if nameprefix == self.default_prefix: 
5679              raise _Error('cannot save to default name "' + nameprefix + '...", chose another name') 
5680   
5681          if nameprefix == self.name_prefix: 
5682              return 
5683   
5684          for name in CMADataLogger.names: 
5685              open(nameprefix + name + '.dat', 'w').write(open(self.name_prefix + name + '.dat').read()) 
5686   
5687          if switch: 
5688              self.name_prefix = nameprefix 
5689   
5690 -    def plot(self, fig=None, iabscissa=1, iteridx=None, 
5691               plot_mean=False, # TODO: plot_mean default should be False 
5692               foffset=1e-19, x_opt=None, fontsize=10): 
5693          """ 
5694          plot data from a `DEAPCMADataLogger` (using the files written by the logger). 
5695   
5696          Arguments 
5697          --------- 
5698              `fig` 
5699                  figure number, by default 325 
5700              `iabscissa` 
5701                  ``0==plot`` versus iteration count, 
5702                  ``1==plot`` versus function evaluation number 
5703              `iteridx` 
5704                  iteration indices to plot 
5705   
5706          Return `CMADataLogger` itself. 
5707   
5708          Examples 
5709          -------- 
5710          :: 
5711   
5712              import cma 
5713              logger = cma.CMADataLogger()  # with default name 
5714              # try to plot the "default logging" data (e.g. from previous fmin calls) 
5715              logger.plot() # to continue you might need to close the pop-up window 
5716                            # once and call plot() again. 
5717                            # This behavior seems to disappear in subsequent 
5718                            # calls of plot(). Also using ipython with -pylab 
5719                            # option might help. 
5720              cma.savefig('fig325.png')  # save current figure 
5721              logger.closefig() 
5722   
5723          Dependencies: matlabplotlib/pyplot. 
5724   
5725          """ 
5726   
5727          dat = self.load(self.name_prefix) 
5728   
5729          try: 
5730              # pyplot: prodedural interface for matplotlib 
5731              from  matplotlib.pyplot import figure, ioff, ion, subplot, semilogy, hold, plot, grid, \ 
5732                   axis, title, text, xlabel, isinteractive, draw, gcf 
5733   
5734          except ImportError: 
5735              ImportError('could not find matplotlib.pyplot module, function plot() is not available') 
5736              return 
5737   
5738          if fontsize and pyplot.rcParams['font.size'] != fontsize: 
5739              print('global variable pyplot.rcParams[\'font.size\'] set (from ' + 
5740                    str(pyplot.rcParams['font.size']) + ') to ' + str(fontsize)) 
5741              pyplot.rcParams['font.size'] = fontsize  # subtracted in the end, but return can happen inbetween 
5742   
5743          if fig: 
5744              figure(fig) 
5745          else: 
5746              figure(325) 
5747              # show()  # should not be necessary 
5748          self.fighandle = gcf()  # fighandle.number 
5749   
5750          if iabscissa not in (0, 1): 
5751              iabscissa = 1 
5752          interactive_status = isinteractive() 
5753          ioff()  # prevents immediate drawing 
5754   
5755          dat.x = dat.xmean  # this is the genotyp 
5756          if not plot_mean: 
5757              if len(dat.x) < 2: 
5758                  print('not enough data to plot recent x, using mean instead') 
5759              else: 
5760                  dat.x = dat.xrecent 
5761          if iteridx is not None: 
5762              dat.f = dat.f[np.where([x in iteridx for x in dat.f[:, 0]])[0], :] 
5763              dat.D = dat.D[np.where([x in iteridx for x in dat.D[:, 0]])[0], :] 
5764              iteridx.append(dat.x[-1, 1])  # last entry is artificial 
5765              dat.x = dat.x[np.where([x in iteridx for x in dat.x[:, 0]])[0], :] 
5766              dat.std = dat.std[np.where([x in iteridx for x in dat.std[:, 0]])[0], :] 
5767   
5768          if iabscissa == 0: 
5769              xlab = 'iterations' 
5770          elif iabscissa == 1: 
5771              xlab = 'function evaluations' 
5772   
5773          # use fake last entry in x and std for line extension-annotation 
5774          if dat.x.shape[1] < 100: 
5775              minxend = int(1.06 * dat.x[-2, iabscissa]) 
5776              # write y-values for individual annotation into dat.x 
5777              dat.x[-1, iabscissa] = minxend  # TODO: should be ax[1] 
5778              idx = np.argsort(dat.x[-2, 5:]) 
5779              idx2 = np.argsort(idx) 
5780              if x_opt is None: 
5781                  dat.x[-1, 5 + idx] = np.linspace(np.min(dat.x[:, 5:]), 
5782                              np.max(dat.x[:, 5:]), dat.x.shape[1] - 5) 
5783              else: 
5784                  dat.x[-1, 5 + idx] = np.logspace(np.log10(np.min(abs(dat.x[:, 5:]))), 
5785                              np.log10(np.max(abs(dat.x[:, 5:]))), dat.x.shape[1] - 5) 
5786          else: 
5787              minxend = 0 
5788   
5789          if len(dat.f) == 0: 
5790              print('nothing to plot') 
5791              return 
5792   
5793          # not in use anymore, see formatter above 
5794          # xticklocs = np.arange(5) * np.round(minxend/4., -int(np.log10(minxend/4.))) 
5795   
5796          # dfit(dfit<1e-98) = NaN; 
5797   
5798          ioff()  # turns update off 
5799   
5800          # TODO: if abscissa==0 plot in chunks, ie loop over subsets where dat.f[:,0]==countiter is monotonous 
5801   
5802          subplot(2, 2, 1) 
5803          self.plotdivers(dat, iabscissa, foffset) 
5804   
5805          # TODO: modularize also the remaining subplots 
5806          subplot(2, 2, 2) 
5807          hold(False) 
5808          if x_opt is not None:  # TODO: differentate neg and pos? 
5809              semilogy(dat.x[:, iabscissa], abs(dat.x[:, 5:]) - x_opt, '-') 
5810          else: 
5811              plot(dat.x[:, iabscissa], dat.x[:, 5:], '-') 
5812          hold(True) 
5813          grid(True) 
5814          ax = array(axis()) 
5815          # ax[1] = max(minxend, ax[1]) 
5816          axis(ax) 
5817          ax[1] -= 1e-6 
5818          if dat.x.shape[1] < 100: 
5819              yy = np.linspace(ax[2] + 1e-6, ax[3] - 1e-6, dat.x.shape[1] - 5) 
5820              # yyl = np.sort(dat.x[-1,5:]) 
5821              idx = np.argsort(dat.x[-1, 5:]) 
5822              idx2 = np.argsort(idx) 
5823              if x_opt is not None: 
5824                  semilogy([dat.x[-1, iabscissa], ax[1]], [abs(dat.x[-1, 5:]), yy[idx2]], 'k-')  # line from last data point 
5825                  semilogy(np.dot(dat.x[-2, iabscissa], [1, 1]), array([ax[2] + 1e-6, ax[3] - 1e-6]), 'k-') 
5826              else: 
5827                  # plot([dat.x[-1, iabscissa], ax[1]], [dat.x[-1,5:], yy[idx2]], 'k-') # line from last data point 
5828                  plot(np.dot(dat.x[-2, iabscissa], [1, 1]), array([ax[2] + 1e-6, ax[3] - 1e-6]), 'k-') 
5829              # plot(array([dat.x[-1, iabscissa], ax[1]]), 
5830              #      reshape(array([dat.x[-1,5:], yy[idx2]]).flatten(), (2,4)), '-k') 
5831              for i in range(len(idx)): 
5832                  # TODOqqq: annotate phenotypic value!? 
5833                  # text(ax[1], yy[i], 'x(' + str(idx[i]) + ')=' + str(dat.x[-2,5+idx[i]])) 
5834                  text(dat.x[-1, iabscissa], dat.x[-1, 5 + i], 'x(' + str(i) + ')=' + str(dat.x[-2, 5 + i])) 
5835   
5836          i = 2  # find smallest i where iteration count differs (in case the same row appears twice) 
5837          while i < len(dat.f) and dat.f[-i][0] == dat.f[-1][0]: 
5838              i += 1 
5839          title('Object Variables (' + ('mean' if plot_mean else 'curr best') + 
5840                  ', ' + str(dat.x.shape[1] - 5) + '-D, popsize~' + 
5841                  (str(int((dat.f[-1][1] - dat.f[-i][1]) / (dat.f[-1][0] - dat.f[-i][0]))) 
5842                      if len(dat.f.T[0]) > 1 and dat.f[-1][0] > dat.f[-i][0] else 'NA') 
5843                  + ')') 
5844          # pyplot.xticks(xticklocs) 
5845   
5846          # Scaling 
5847          subplot(2, 2, 3) 
5848          hold(False) 
5849          semilogy(dat.D[:, iabscissa], dat.D[:, 5:], '-b') 
5850          hold(True) 
5851          grid(True) 
5852          ax = array(axis()) 
5853          # ax[1] = max(minxend, ax[1]) 
5854          axis(ax) 
5855          title('Scaling (All Main Axes)') 
5856          # pyplot.xticks(xticklocs) 
5857          xlabel(xlab) 
5858   
5859          # standard deviations 
5860          subplot(2, 2, 4) 
5861          hold(False) 
5862          # remove sigma from stds (graphs become much better readible) 
5863          dat.std[:, 5:] = np.transpose(dat.std[:, 5:].T / dat.std[:, 2].T) 
5864          # ax = array(axis()) 
5865          # ax[1] = max(minxend, ax[1]) 
5866          # axis(ax) 
5867          if 1 < 2 and dat.std.shape[1] < 100: 
5868              # use fake last entry in x and std for line extension-annotation 
5869              minxend = int(1.06 * dat.x[-2, iabscissa]) 
5870              dat.std[-1, iabscissa] = minxend  # TODO: should be ax[1] 
5871              idx = np.argsort(dat.std[-2, 5:]) 
5872              idx2 = np.argsort(idx) 
5873              dat.std[-1, 5 + idx] = np.logspace(np.log10(np.min(dat.std[:, 5:])), 
5874                              np.log10(np.max(dat.std[:, 5:])), dat.std.shape[1] - 5) 
5875   
5876              dat.std[-1, iabscissa] = minxend  # TODO: should be ax[1] 
5877              yy = np.logspace(np.log10(ax[2]), np.log10(ax[3]), dat.std.shape[1] - 5) 
5878              # yyl = np.sort(dat.std[-1,5:]) 
5879              idx = np.argsort(dat.std[-1, 5:]) 
5880              idx2 = np.argsort(idx) 
5881              # plot(np.dot(dat.std[-2, iabscissa],[1,1]), array([ax[2]+1e-6, ax[3]-1e-6]), 'k-') # vertical separator 
5882              # vertical separator 
5883              plot(np.dot(dat.std[-2, iabscissa], [1, 1]), array([np.min(dat.std[-2, 5:]), np.max(dat.std[-2, 5:])]), 'k-') 
5884              hold(True) 
5885              # plot([dat.std[-1, iabscissa], ax[1]], [dat.std[-1,5:], yy[idx2]], 'k-') # line from last data point 
5886              for i in rglen((idx)): 
5887                  # text(ax[1], yy[i], ' '+str(idx[i])) 
5888                  text(dat.std[-1, iabscissa], dat.std[-1, 5 + i], ' ' + str(i)) 
5889          semilogy(dat.std[:, iabscissa], dat.std[:, 5:], '-') 
5890          grid(True) 
5891          title('Standard Deviations in All Coordinates') 
5892          # pyplot.xticks(xticklocs) 
5893          xlabel(xlab) 
5894          ion() 
5895          draw()  # does not suffice 
5896          show() 
5897   
5898          return self 
5899   
5900   
5901      # ____________________________________________________________ 
5902      # ____________________________________________________________ 
5903      # 
5904      @staticmethod 
5905 -    def plotdivers(dat, iabscissa, foffset): 
5906          """helper function for `plot()` that plots all what is 
5907          in the upper left subplot like fitness, sigma, etc. 
5908   
5909          Arguments 
5910          --------- 
5911              `iabscissa` in ``(0,1)`` 
5912                  0==versus fevals, 1==versus iteration 
5913              `foffset` 
5914                  offset to fitness for log-plot 
5915   
5916           :See: `plot()` 
5917   
5918          """ 
5919          from  matplotlib.pyplot import semilogy, hold, grid, \ 
5920                   axis, title, text 
5921          fontsize = pyplot.rcParams['font.size'] 
5922   
5923          hold(False) 
5924   
5925          dfit = dat.f[:, 5] - min(dat.f[:, 5]) 
5926          dfit[dfit < 1e-98] = np.NaN 
5927   
5928          if dat.f.shape[1] > 7: 
5929              # semilogy(dat.f[:, iabscissa], abs(dat.f[:,[6, 7, 10, 12]])+foffset,'-k') 
5930              semilogy(dat.f[:, iabscissa], abs(dat.f[:, [6, 7]]) + foffset, '-k') 
5931              hold(True) 
5932   
5933          # (larger indices): additional fitness data, for example constraints values 
5934          if dat.f.shape[1] > 8: 
5935              # dd = abs(dat.f[:,7:]) + 10*foffset 
5936              # dd = np.where(dat.f[:,7:]==0, np.NaN, dd) # cannot be 
5937              semilogy(dat.f[:, iabscissa], np.abs(dat.f[:, 8:]) + 10 * foffset, 'm') 
5938              hold(True) 
5939   
5940          idx = np.where(dat.f[:, 5] > 1e-98)[0]  # positive values 
5941          semilogy(dat.f[idx, iabscissa], dat.f[idx, 5] + foffset, '.b') 
5942          hold(True) 
5943          grid(True) 
5944   
5945          idx = np.where(dat.f[:, 5] < -1e-98)  # negative values 
5946          semilogy(dat.f[idx, iabscissa], abs(dat.f[idx, 5]) + foffset, '.r') 
5947   
5948          semilogy(dat.f[:, iabscissa], abs(dat.f[:, 5]) + foffset, '-b') 
5949          semilogy(dat.f[:, iabscissa], dfit, '-c') 
5950   
5951          # overall minimum 
5952          i = np.argmin(dat.f[:, 5]) 
5953          semilogy(dat.f[i, iabscissa] * np.ones(2), dat.f[i, 5] * np.ones(2), 'rd') 
5954          # semilogy(dat.f[-1, iabscissa]*np.ones(2), dat.f[-1,4]*np.ones(2), 'rd') 
5955   
5956          # AR and sigma 
5957          semilogy(dat.f[:, iabscissa], dat.f[:, 3], '-r')  # AR 
5958          semilogy(dat.f[:, iabscissa], dat.f[:, 2], '-g')  # sigma 
5959          semilogy(dat.std[:-1, iabscissa], np.vstack([list(map(max, dat.std[:-1, 5:])), list(map(min, dat.std[:-1, 5:]))]).T, 
5960                       '-m', linewidth=2) 
5961          text(dat.std[-2, iabscissa], max(dat.std[-2, 5:]), 'max std', fontsize=fontsize) 
5962          text(dat.std[-2, iabscissa], min(dat.std[-2, 5:]), 'min std', fontsize=fontsize) 
5963          ax = array(axis()) 
5964          # ax[1] = max(minxend, ax[1]) 
5965          axis(ax) 
5966          text(ax[0] + 0.01, ax[2],  # 10**(log10(ax[2])+0.05*(log10(ax[3])-log10(ax[2]))), 
5967               '.f_recent=' + repr(dat.f[-1, 5])) 
5968   
5969          # title('abs(f) (blue), f-min(f) (cyan), Sigma (green), Axis Ratio (red)') 
5970          title('blue:abs(f), cyan:f-min(f), green:sigma, red:axis ratio', fontsize=fontsize - 1) 
5971          # pyplot.xticks(xticklocs) 
5972   
5973   
5974   
5975 -    def downsampling(self, factor=10, first=3, switch=True): 
5976          """ 
5977          rude downsampling of a `CMADataLogger` data file by `factor`, keeping 
5978          also the first `first` entries. This function is a stump and subject 
5979          to future changes. 
5980   
5981          Arguments 
5982          --------- 
5983             - `factor` -- downsampling factor 
5984             - `first` -- keep first `first` entries 
5985             - `switch` -- switch the new logger name to oldname+'down' 
5986   
5987          Details 
5988          ------- 
5989          ``self.name_prefix+'down'`` files are written 
5990   
5991          Example 
5992          ------- 
5993          :: 
5994   
5995              import cma 
5996              cma.downsampling()  # takes outcmaes* files 
5997              cma.plot('outcmaesdown') 
5998   
5999          """ 
6000          newprefix = self.name_prefix + 'down' 
6001          for name in CMADataLogger.names: 
6002              f = open(newprefix + name + '.dat', 'w') 
6003              iline = 0 
6004              cwritten = 0 
6005              for line in open(self.name_prefix + name + '.dat'): 
6006                  if iline < first or iline % factor == 0: 
6007                      f.write(line) 
6008                      cwritten += 1 
6009                  iline += 1 
6010              f.close() 
6011              print('%d' % (cwritten) + ' lines written in ' + newprefix + name + '.dat') 
6012          if switch: 
6013              self.name_prefix += 'down' 
6014          return self 
6015   
6016      # ____________________________________________________________ 
6017      # ____________________________________________________________ 
6018      # 
6019 -    def disp_header(self): 
6020          heading = 'Iterat Nfevals  function value    axis ratio maxstd   minstd' 
6021          print(heading) 
6022   
6023 -    def disp(self, idx=100):  # r_[0:5,1e2:1e9:1e2,-10:0]): 
6024          """displays selected data from (files written by) the class `CMADataLogger`. 
6025   
6026          Arguments 
6027          --------- 
6028             `idx` 
6029                 indices corresponding to rows in the data file; 
6030                 if idx is a scalar (int), the first two, then every idx-th, 
6031                 and the last three rows are displayed. Too large index values are removed. 
6032                 If ``len(idx) == 1``, only a single row is displayed, e.g. the last 
6033                 entry when ``idx == [-1]``. 
6034   
6035          Example 
6036          ------- 
6037          >>> import cma, numpy as np 
6038          >>> res = cma.fmin(cma.fcts.elli, 7 * [0.1], 1, {'verb_disp':1e9})  # generate data 
6039          >>> assert res[1] < 1e-9 
6040          >>> assert res[2] < 4400 
6041          >>> l = cma.CMADataLogger()  # == res[-1], logger with default name, "points to" above data 
6042          >>> l.disp([0,-1])  # first and last 
6043          >>> l.disp(20)  # some first/last and every 20-th line 
6044          >>> l.disp(np.r_[0:999999:100, -1]) # every 100-th and last 
6045          >>> l.disp(np.r_[0, -10:0]) # first and ten last 
6046          >>> cma.disp(l.name_prefix, np.r_[0::100, -10:])  # the same as l.disp(...) 
6047   
6048          Details 
6049          ------- 
6050          The data line with the best f-value is displayed as last line. 
6051   
6052          :See: `disp()` 
6053   
6054          """ 
6055   
6056          filenameprefix = self.name_prefix 
6057   
6058          def printdatarow(dat, iteration): 
6059              """print data of iteration i""" 
6060              i = np.where(dat.f[:, 0] == iteration)[0][0] 
6061              j = np.where(dat.std[:, 0] == iteration)[0][0] 
6062              print('%5d' % (int(dat.f[i, 0])) + ' %6d' % (int(dat.f[i, 1])) + ' %.14e' % (dat.f[i, 5]) + 
6063                    ' %5.1e' % (dat.f[i, 3]) + 
6064                    ' %6.2e' % (max(dat.std[j, 5:])) + ' %6.2e' % min(dat.std[j, 5:])) 
6065   
6066          dat = CMADataLogger(filenameprefix).load() 
6067          ndata = dat.f.shape[0] 
6068   
6069          # map index to iteration number, is difficult if not all iteration numbers exist 
6070          # idx = idx[np.where(map(lambda x: x in dat.f[:,0], idx))[0]] # TODO: takes pretty long 
6071          # otherwise: 
6072          if idx is None: 
6073              idx = 100 
6074          if np.isscalar(idx): 
6075              # idx = np.arange(0, ndata, idx) 
6076              if idx: 
6077                  idx = np.r_[0, 1, idx:ndata - 3:idx, -3:0] 
6078              else: 
6079                  idx = np.r_[0, 1, -3:0] 
6080   
6081          idx = array(idx) 
6082          idx = idx[idx <= ndata]  # TODO: shouldn't this be "<"? 
6083          idx = idx[-idx <= ndata] 
6084          iters = dat.f[idx, 0] 
6085          idxbest = np.argmin(dat.f[:, 5]) 
6086          iterbest = dat.f[idxbest, 0] 
6087          if len(iters) == 1: 
6088              printdatarow(dat, iters[0]) 
6089          else: 
6090              self.disp_header() 
6091              for i in iters: 
6092                  printdatarow(dat, i) 
6093              self.disp_header() 
6094              printdatarow(dat, iterbest) 
6095          sys.stdout.flush() 
6096   
6097  last_figure_number = 324 
6098 -def plot(name=None, fig=None, abscissa=1, iteridx=None, 
6099           plot_mean=False, 
6100           foffset=1e-19, x_opt=None, fontsize=10): 
6101      """ 
6102      plot data from files written by a `CMADataLogger`, 
6103      the call ``cma.plot(name, **argsdict)`` is a shortcut for 
6104      ``cma.CMADataLogger(name).plot(**argsdict)`` 
6105   
6106      Arguments 
6107      --------- 
6108          `name` 
6109              name of the logger, filename prefix, None evaluates to 
6110              the default 'outcmaes' 
6111          `fig` 
6112              filename or figure number, or both as a tuple (any order) 
6113          `abscissa` 
6114              0==plot versus iteration count, 
6115              1==plot versus function evaluation number 
6116          `iteridx` 
6117              iteration indices to plot 
6118   
6119      Return `None` 
6120   
6121      Examples 
6122      -------- 
6123      :: 
6124   
6125         cma.plot();  # the optimization might be still 
6126                      # running in a different shell 
6127         cma.savefig('fig325.png') 
6128         cma.closefig() 
6129   
6130         cdl = cma.CMADataLogger().downsampling().plot() 
6131         # in case the file sizes are large 
6132   
6133      Details 
6134      ------- 
6135      Data from codes in other languages (C, Java, Matlab, Scilab) have the same 
6136      format and can be plotted just the same. 
6137   
6138      :See: `CMADataLogger`, `CMADataLogger.plot()` 
6139   
6140      """ 
6141      global last_figure_number 
6142      if not fig: 
6143          last_figure_number += 1 
6144          fig = last_figure_number 
6145      last_figure_number = fig 
6146      CMADataLogger(name).plot(fig, abscissa, iteridx, plot_mean, foffset, 
6147                               x_opt, fontsize) 
6148   
6149 -def disp(name=None, idx=None): 
6150      """displays selected data from (files written by) the class `CMADataLogger`. 
6151   
6152      The call ``cma.disp(name, idx)`` is a shortcut for ``cma.CMADataLogger(name).disp(idx)``. 
6153   
6154      Arguments 
6155      --------- 
6156          `name` 
6157              name of the logger, filename prefix, `None` evaluates to 
6158              the default ``'outcmaes'`` 
6159          `idx` 
6160              indices corresponding to rows in the data file; by 
6161              default the first five, then every 100-th, and the last 
6162              10 rows. Too large index values are removed. 
6163   
6164      Examples 
6165      -------- 
6166      :: 
6167   
6168         import cma, numpy 
6169         # assume some data are available from previous runs 
6170         cma.disp(None,numpy.r_[0,-1])  # first and last 
6171         cma.disp(None,numpy.r_[0:1e9:100,-1]) # every 100-th and last 
6172         cma.disp(idx=numpy.r_[0,-10:0]) # first and ten last 
6173         cma.disp(idx=numpy.r_[0:1e9:1e3,-10:0]) 
6174   
6175      :See: `CMADataLogger.disp()` 
6176   
6177      """ 
6178      return CMADataLogger(name if name else 'outcmaes' 
6179                           ).disp(idx) 
6180   
6181  # ____________________________________________________________ 
6182 -def _fileToMatrix(file_name): 
6183      """rudimentary method to read in data from a file""" 
6184      # TODO: np.loadtxt() might be an alternative 
6185      #     try: 
6186      if 1 < 3: 
6187          lres = [] 
6188          for line in open(file_name, 'r').readlines(): 
6189              if len(line) > 0 and line[0] not in ('%', '#'): 
6190                  lres.append(list(map(float, line.split()))) 
6191          res = lres 
6192      while res != [] and res[0] == []:  # remove further leading empty lines 
6193          del res[0] 
6194      return res 
6195      #     except: 
6196      print('could not read file ' + file_name) 
6197   
6198  # ____________________________________________________________ 
6199  # ____________________________________________________________ 
6200 -class NoiseHandler(object): 
6201      """Noise handling according to [Hansen et al 2009, A Method for 
6202      Handling Uncertainty in Evolutionary Optimization...] 
6203   
6204      The interface of this class is yet versatile and subject to changes. 
6205   
6206      The attribute ``evaluations`` serves to control the noise via 
6207      number of evaluations, for example in ``fmin`` or with 
6208      `ask_and_eval()`. The parameter ``maxevals`` (second parameter) 
6209      provides the upper bound, or lower and upper bound, or lower and 
6210      upper bound and initial value, all 1 by default, compare also the 
6211      second example. 
6212   
6213      Examples 
6214      -------- 
6215      Minimal example together with `fmin` on a non-noisy function: 
6216   
6217      >>> import cma 
6218      >>> cma.fmin(cma.felli, 7 * [1], 1, noise_handler=cma.NoiseHandler(7))  # dimension 7 
6219   
6220      More verbose example in the optimization loop with a noisy function 
6221      defined in ``func``: 
6222   
6223      >>> import cma, numpy as np 
6224      >>> func = lambda x: cma.fcts.sphere(x) * (1 + 4 * np.random.randn() / len(x))  # cma.Fcts.noisysphere 
6225      >>> es = cma.CMAEvolutionStrategy(np.ones(10), 1) 
6226      >>> nh = cma.NoiseHandler(es.N, maxevals=[1, 1, 30]) 
6227      >>> while not es.stop(): 
6228      ...     X, fit_vals = es.ask_and_eval(func, evaluations=nh.evaluations) 
6229      ...     es.tell(X, fit_vals)  # prepare for next iteration 
6230      ...     es.sigma *= nh(X, fit_vals, func, es.ask)  # see method __call__ 
6231      ...     es.countevals += nh.evaluations_just_done  # this is a hack, not important though 
6232      ...     es.logger.add(more_data = [nh.evaluations, nh.noiseS])  # add a data point 
6233      ...     es.disp() 
6234      ...     # nh.maxevals = ...  it might be useful to start with smaller values and then increase 
6235      >>> print(es.stop()) 
6236      >>> print(es.result()[-2])  # take mean value, the best solution is totally off 
6237      >>> assert sum(es.result()[-2]**2) < 1e-9 
6238      >>> print(X[np.argmin(fit_vals)])  # not bad, but probably worse than the mean 
6239      >>> # es.logger.plot() 
6240   
6241   
6242      The command ``logger.plot()`` will plot the logged data. 
6243   
6244      The noise options of `fmin()` control a `NoiseHandler` instance similar to this 
6245      example. The command ``cma.CMAOptions('noise')`` lists in effect the parameters of 
6246      `__init__` apart from ``aggregate``. 
6247   
6248      Details 
6249      ------- 
6250      The parameters reevals, theta, c_s, and alpha_t are set differently 
6251      than in the original publication, see method `__init__()`. For a 
6252      very small population size, say popsize <= 5, the measurement 
6253      technique based on rank changes is likely to fail. 
6254   
6255      Missing Features 
6256      ---------------- 
6257      In case no noise is found, ``self.lam_reeval`` should be adaptive 
6258      and get at least as low as 1 (however the possible savings from this 
6259      are rather limited). Another option might be to decide during the 
6260      first call by a quantitative analysis of fitness values whether 
6261      ``lam_reeval`` is set to zero. More generally, an automatic noise 
6262      mode detection might also set the covariance matrix learning rates 
6263      to smaller values. 
6264   
6265      :See: `fmin()`, `ask_and_eval()` 
6266   
6267      """ 
6268      # TODO: for const additive noise a better version might be with alphasigma also used for sigma-increment, 
6269      # while all other variance changing sources are removed (because they are intrinsically biased). Then 
6270      # using kappa to get convergence (with unit sphere samples): noiseS=0 leads to a certain kappa increasing rate? 
6271 -    def __init__(self, N, maxevals=[1, 1, 1], aggregate=np.median, reevals=None, epsilon=1e-7, parallel=False): 
6272          """parameters are 
6273   
6274              `N` 
6275                  dimension, (only) necessary to adjust the internal "alpha"-parameters 
6276              `maxevals` 
6277                  maximal value for ``self.evaluations``, where 
6278                  ``self.evaluations`` function calls are aggregated for 
6279                  noise treatment. With ``maxevals == 0`` the noise 
6280                  handler is (temporarily) "switched off". If `maxevals` 
6281                  is a list, min value and (for >2 elements) median are 
6282                  used to define minimal and initial value of 
6283                  ``self.evaluations``. Choosing ``maxevals > 1`` is only 
6284                  reasonable, if also the original ``fit`` values (that 
6285                  are passed to `__call__`) are computed by aggregation of 
6286                  ``self.evaluations`` values (otherwise the values are 
6287                  not comparable), as it is done within `fmin()`. 
6288              `aggregate` 
6289                  function to aggregate single f-values to a 'fitness', e.g. 
6290                  ``np.median``. 
6291              `reevals` 
6292                  number of solutions to be reevaluated for noise measurement, 
6293                  can be a float, by default set to ``2 + popsize/20``, where 
6294                  ``popsize = len(fit)`` in ``__call__``. 
6295                  zero switches noise handling off. 
6296              `epsilon` 
6297                  multiplier for perturbation of the reevaluated solutions 
6298              `parallel` 
6299                  a single f-call with all resampled solutions 
6300   
6301              :See: `fmin()`, `CMAOptions`, `CMAEvolutionStrategy.ask_and_eval()` 
6302   
6303          """ 
6304          self.lam_reeval = reevals  # 2 + popsize/20, see method indices(), originally 2 + popsize/10 
6305          self.epsilon = epsilon 
6306          self.parallel = parallel 
6307          self.theta = 0.5  # originally 0.2 
6308          self.cum = 0.3  # originally 1, 0.3 allows one disagreement of current point with resulting noiseS 
6309          self.alphasigma = 1 + 2 / (N + 10) # unit sphere sampling: 1 + 1 / (N + 10) 
6310          self.alphaevals = 1 + 2 / (N + 10)  # originally 1.5 
6311          self.alphaevalsdown = self.alphaevals**-0.25  # originally 1/1.5 
6312          # zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz 
6313          self.evaluations = 1  # to aggregate for a single f-evaluation 
6314          self.minevals = 1 
6315          self.maxevals = int(np.max(maxevals)) 
6316          if hasattr(maxevals, '__contains__'):  # i.e. can deal with ``in`` 
6317              if len(maxevals) > 1: 
6318                  self.minevals = min(maxevals) 
6319                  self.evaluations = self.minevals 
6320              if len(maxevals) > 2: 
6321                  self.evaluations = np.median(maxevals) 
6322          self.f_aggregate = aggregate 
6323          self.evaluations_just_done = 0  # actually conducted evals, only for documentation 
6324          self.noiseS = 0 
6325   
6326 -    def __call__(self, X, fit, func, ask=None, args=()): 
6327          """proceed with noise measurement, set anew attributes ``evaluations`` 
6328          (proposed number of evaluations to "treat" noise) and ``evaluations_just_done`` 
6329          and return a factor for increasing sigma. 
6330   
6331          Parameters 
6332          ---------- 
6333              `X` 
6334                  a list/sequence/vector of solutions 
6335              `fit` 
6336                  the respective list of function values 
6337              `func` 
6338                  the objective function, ``fit[i]`` corresponds to ``func(X[i], *args)`` 
6339              `ask` 
6340                  a method to generate a new, slightly disturbed solution. The argument 
6341                  is (only) mandatory if ``epsilon`` is not zero, see `__init__()`. 
6342              `args` 
6343                  optional additional arguments to `func` 
6344   
6345          Details 
6346          ------- 
6347          Calls the methods ``reeval()``, ``update_measure()`` and ``treat()`` in this order. 
6348          ``self.evaluations`` is adapted within the method `treat()`. 
6349   
6350          """ 
6351          self.evaluations_just_done = 0 
6352          if not self.maxevals or self.lam_reeval == 0: 
6353              return 1.0 
6354          res = self.reeval(X, fit, func, ask, args) 
6355          if not len(res): 
6356              return 1.0 
6357          self.update_measure() 
6358          return self.treat() 
6359   
6360 -    def get_evaluations(self): 
6361          """return ``self.evaluations``, the number of evalutions to get a single fitness measurement""" 
6362          return self.evaluations 
6363   
6364 -    def treat(self): 
6365          """adapt self.evaluations depending on the current measurement value 
6366          and return ``sigma_fac in (1.0, self.alphasigma)`` 
6367   
6368          """ 
6369          if self.noiseS > 0: 
6370              self.evaluations = min((self.evaluations * self.alphaevals, self.maxevals)) 
6371              return self.alphasigma 
6372          else: 
6373              self.evaluations = max((self.evaluations * self.alphaevalsdown, self.minevals)) 
6374              return 1.0  # / self.alphasigma 
6375   
6376 -    def reeval(self, X, fit, func, ask, args=()): 
6377          """store two fitness lists, `fit` and ``fitre`` reevaluating some 
6378          solutions in `X`. 
6379          ``self.evaluations`` evaluations are done for each reevaluated 
6380          fitness value. 
6381          See `__call__()`, where `reeval()` is called. 
6382   
6383          """ 
6384          self.fit = list(fit) 
6385          self.fitre = list(fit) 
6386          self.idx = self.indices(fit) 
6387          if not len(self.idx): 
6388              return self.idx 
6389          evals = int(self.evaluations) if self.f_aggregate else 1 
6390          fagg = np.median if self.f_aggregate is None else self.f_aggregate 
6391          for i in self.idx: 
6392              X_i = X[i] 
6393              if self.epsilon: 
6394                  if self.parallel: 
6395                      self.fitre[i] = fagg(func(ask(evals, X_i, self.epsilon), *args)) 
6396                  else: 
6397                      self.fitre[i] = fagg([func(ask(1, X_i, self.epsilon)[0], *args) 
6398                                              for _k in xrange(evals)]) 
6399              else: 
6400                  self.fitre[i] = fagg([func(X_i, *args) for _k in xrange(evals)]) 
6401          self.evaluations_just_done = evals * len(self.idx) 
6402          return self.fit, self.fitre, self.idx 
6403   
6404 -    def update_measure(self): 
6405          """updated noise level measure using two fitness lists ``self.fit`` and 
6406          ``self.fitre``, return ``self.noiseS, all_individual_measures``. 
6407   
6408          Assumes that `self.idx` contains the indices where the fitness 
6409          lists differ 
6410   
6411          """ 
6412          lam = len(self.fit) 
6413          idx = np.argsort(self.fit + self.fitre) 
6414          ranks = np.argsort(idx).reshape((2, lam)) 
6415          rankDelta = ranks[0] - ranks[1] - np.sign(ranks[0] - ranks[1]) 
6416   
6417          # compute rank change limits using both ranks[0] and ranks[1] 
6418          r = np.arange(1, 2 * lam)  # 2 * lam - 2 elements 
6419          limits = [0.5 * (Mh.prctile(np.abs(r - (ranks[0, i] + 1 - (ranks[0, i] > ranks[1, i]))), 
6420                                        self.theta * 50) + 
6421                           Mh.prctile(np.abs(r - (ranks[1, i] + 1 - (ranks[1, i] > ranks[0, i]))), 
6422                                        self.theta * 50)) 
6423                      for i in self.idx] 
6424          # compute measurement 
6425          #                               max: 1 rankchange in 2*lambda is always fine 
6426          s = np.abs(rankDelta[self.idx]) - Mh.amax(limits, 1)  # lives roughly in 0..2*lambda 
6427          self.noiseS += self.cum * (np.mean(s) - self.noiseS) 
6428          return self.noiseS, s 
6429   
6430 -    def indices(self, fit): 
6431          """return the set of indices to be reevaluated for noise measurement, 
6432          taking the ``lam_reeval`` best from the first ``2 * lam_reeval + 2`` 
6433          values. 
6434   
6435          Given the first values are the earliest, this is a useful policy also 
6436          with a time changing objective. 
6437   
6438          """ 
6439          lam = self.lam_reeval if self.lam_reeval else 2 + len(fit) / 20 
6440          reev = int(lam) + ((lam % 1) > np.random.rand()) 
6441          return np.argsort(array(fit, copy=False)[:2 * (reev + 1)])[:reev] 
6442   
6443  # ____________________________________________________________ 
6444  # ____________________________________________________________ 
6445 -class Sections(object): 
6446      """plot sections through an objective function. 
6447   
6448      A first rational thing to do, when facing an (expensive) 
6449      application. By default 6 points in each coordinate are evaluated. 
6450      This class is still experimental. 
6451   
6452      Examples 
6453      -------- 
6454   
6455      >>> import cma, numpy as np 
6456      >>> s = cma.Sections(cma.Fcts.rosen, np.zeros(3)).do(plot=False) 
6457      >>> s.do(plot=False)  # evaluate the same points again, i.e. check for noise 
6458      >> try: 
6459      ...     s.plot() 
6460      ... except: 
6461      ...     print('plotting failed: matplotlib.pyplot package is missing?') 
6462   
6463      Details 
6464      ------- 
6465      Data are saved after each function call during `do()`. The filename 
6466      is attribute ``name`` and by default ``str(func)``, see `__init__()`. 
6467   
6468      A random (orthogonal) basis can be generated with 
6469      ``cma.Rotation()(np.eye(3))``. 
6470   
6471      CAVEAT: The default name is unique in the function name, but it 
6472      should be unique in all parameters of `__init__()` but `plot_cmd` 
6473      and `load`. If, for example, a different basis is chosen, either 
6474      the name must be changed or the ``.pkl`` file containing the 
6475      previous data must first be renamed or deleted. 
6476   
6477      ``s.res`` is a dictionary with an entry for each "coordinate" ``i`` 
6478      and with an entry ``'x'``, the middle point. Each entry ``i`` is 
6479      again a dictionary with keys being different dx values and the 
6480      value being a sequence of f-values. For example ``s.res[2][0.1] == 
6481      [0.01, 0.01]``, which is generated using the difference vector ``s 
6482      .basis[2]`` like 
6483   
6484      ``s.res[2][dx] += func(s.res['x'] + dx * s.basis[2])``. 
6485   
6486      :See: `__init__()` 
6487   
6488      """ 
6489 -    def __init__(self, func, x, args=(), basis=None, name=None, 
6490                   plot_cmd=pyplot.plot if pyplot else None, load=True): 
6491          """ 
6492          Parameters 
6493          ---------- 
6494              `func` 
6495                  objective function 
6496              `x` 
6497                  point in search space, middle point of the sections 
6498              `args` 
6499                  arguments passed to `func` 
6500              `basis` 
6501                  evaluated points are ``func(x + locations[j] * basis[i]) for i in len(basis) for j in len(locations)``, 
6502                  see `do()` 
6503              `name` 
6504                  filename where to save the result 
6505              `plot_cmd` 
6506                  command used to plot the data, typically matplotlib pyplots `plot` or `semilogy` 
6507              `load` 
6508                  load previous data from file ``str(func) + '.pkl'`` 
6509   
6510          """ 
6511          self.func = func 
6512          self.args = args 
6513          self.x = x 
6514          self.name = name if name else str(func).replace(' ', '_').replace('>', '').replace('<', '') 
6515          self.plot_cmd = plot_cmd  # or semilogy 
6516          self.basis = np.eye(len(x)) if basis is None else basis 
6517   
6518          try: 
6519              self.load() 
6520              if any(self.res['x'] != x): 
6521                  self.res = {} 
6522                  self.res['x'] = x  # TODO: res['x'] does not look perfect 
6523              else: 
6524                  print(self.name + ' loaded') 
6525          except: 
6526              self.res = {} 
6527              self.res['x'] = x 
6528   
6529 -    def do(self, repetitions=1, locations=np.arange(-0.5, 0.6, 0.2), plot=True): 
6530          """generates, plots and saves function values ``func(y)``, 
6531          where ``y`` is 'close' to `x` (see `__init__()`). The data are stored in 
6532          the ``res`` attribute and the class instance is saved in a file 
6533          with (the weired) name ``str(func)``. 
6534   
6535          Parameters 
6536          ---------- 
6537              `repetitions` 
6538                  for each point, only for noisy functions is >1 useful. For 
6539                  ``repetitions==0`` only already generated data are plotted. 
6540              `locations` 
6541                  coordinated wise deviations from the middle point given in `__init__` 
6542   
6543          """ 
6544          if not repetitions: 
6545              self.plot() 
6546              return 
6547   
6548          res = self.res 
6549          for i in range(len(self.basis)):  # i-th coordinate 
6550              if i not in res: 
6551                  res[i] = {} 
6552              # xx = np.array(self.x) 
6553              # TODO: store res[i]['dx'] = self.basis[i] here? 
6554              for dx in locations: 
6555                  xx = self.x + dx * self.basis[i] 
6556                  xkey = dx  # xx[i] if (self.basis == np.eye(len(self.basis))).all() else dx 
6557                  if xkey not in res[i]: 
6558                      res[i][xkey] = [] 
6559                  n = repetitions 
6560                  while n > 0: 
6561                      n -= 1 
6562                      res[i][xkey].append(self.func(xx, *self.args)) 
6563                      if plot: 
6564                          self.plot() 
6565                      self.save() 
6566          return self 
6567   
6568 -    def plot(self, plot_cmd=None, tf=lambda y: y): 
6569          """plot the data we have, return ``self``""" 
6570          if not plot_cmd: 
6571              plot_cmd = self.plot_cmd 
6572          colors = 'bgrcmyk' 
6573          pyplot.hold(False) 
6574          res = self.res 
6575   
6576          flatx, flatf = self.flattened() 
6577          minf = np.inf 
6578          for i in flatf: 
6579              minf = min((minf, min(flatf[i]))) 
6580          addf = 1e-9 - minf if minf <= 1e-9 else 0 
6581          for i in sorted(res.keys()):  # we plot not all values here 
6582              if isinstance(i, int): 
6583                  color = colors[i % len(colors)] 
6584                  arx = sorted(res[i].keys()) 
6585                  plot_cmd(arx, [tf(np.median(res[i][x]) + addf) for x in arx], color + '-') 
6586                  pyplot.text(arx[-1], tf(np.median(res[i][arx[-1]])), i) 
6587                  pyplot.hold(True) 
6588                  plot_cmd(flatx[i], tf(np.array(flatf[i]) + addf), color + 'o') 
6589          pyplot.ylabel('f + ' + str(addf)) 
6590          pyplot.draw() 
6591          show() 
6592          # raw_input('press return') 
6593          return self 
6594   
6595 -    def flattened(self): 
6596          """return flattened data ``(x, f)`` such that for the sweep through 
6597          coordinate ``i`` we have for data point ``j`` that ``f[i][j] == func(x[i][j])`` 
6598   
6599          """ 
6600          flatx = {} 
6601          flatf = {} 
6602          for i in self.res: 
6603              if isinstance(i, int): 
6604                  flatx[i] = [] 
6605                  flatf[i] = [] 
6606                  for x in sorted(self.res[i]): 
6607                      for d in sorted(self.res[i][x]): 
6608                          flatx[i].append(x) 
6609                          flatf[i].append(d) 
6610          return flatx, flatf 
6611   
6612 -    def save(self, name=None): 
6613          """save to file""" 
6614          import pickle 
6615          name = name if name else self.name 
6616          fun = self.func 
6617          del self.func  # instance method produces error 
6618          pickle.dump(self, open(name + '.pkl', "wb")) 
6619          self.func = fun 
6620          return self 
6621   
6622 -    def load(self, name=None): 
6623          """load from file""" 
6624          import pickle 
6625          name = name if name else self.name 
6626          s = pickle.load(open(name + '.pkl', 'rb')) 
6627          self.res = s.res  # disregard the class 
6628          return self 
6629   
6630  #____________________________________________________________ 
6631  #____________________________________________________________ 
6632 -class _Error(Exception): 
6633      """generic exception of cma module""" 
6634      pass 
6635   
6636  # ____________________________________________________________ 
6637  # ____________________________________________________________ 
6638  # 
6639 -class ElapsedTime(object): 
6640      """using ``time.clock`` with overflow handling to measure CPU time. 
6641   
6642      Example: 
6643   
6644      >>> clock = ElapsedTime()  # clock starts here 
6645      >>> t1 = clock()  # get elapsed CPU time 
6646   
6647      Details: 32-bit C overflows after int(2**32/1e6) == 4294s about 72 min 
6648   
6649      """ 
6650 -    def __init__(self): 
6651          self.tic0 = time.clock() 
6652          self.tic = self.tic0 
6653          self.lasttoc = time.clock() 
6654          self.lastdiff = time.clock() - self.lasttoc 
6655          self.time_to_add = 0 
6656          self.messages = 0 
6657      reset = __init__ 
6658 -    def __call__(self): 
6659          toc = time.clock() 
6660          if toc - self.tic >= self.lasttoc - self.tic: 
6661              self.lastdiff = toc - self.lasttoc 
6662              self.lasttoc = toc 
6663          else:  # overflow, reset self.tic 
6664              if self.messages < 3: 
6665                  self.messages += 1 
6666                  print('  in cma.ElapsedTime: time measure overflow, last difference estimated from', 
6667                          self.tic0, self.tic, self.lasttoc, toc, toc - self.lasttoc, self.lastdiff) 
6668   
6669              self.time_to_add += self.lastdiff + self.lasttoc - self.tic 
6670              self.tic = toc  # reset 
6671              self.lasttoc = toc 
6672          self.elapsedtime = toc - self.tic + self.time_to_add 
6673          return self.elapsedtime 
6674   
6675 -class Misc(object): 
6676      # ____________________________________________________________ 
6677      # ____________________________________________________________ 
6678      # 
6679 -    class MathHelperFunctions(object): 
6680          """static convenience math helper functions, if the function name 
6681          is preceded with an "a", a numpy array is returned 
6682   
6683          """ 
6684          @staticmethod 
6685 -        def aclamp(x, upper): 
6686              return -Misc.MathHelperFunctions.apos(-x, -upper) 
6687          @staticmethod 
6688 -        def equals_approximately(a, b, eps=1e-12): 
6689              if a < 0: 
6690                  a, b = -1 * a, -1 * b 
6691              return (a - eps < b < a + eps) or ((1 - eps) * a < b < (1 + eps) * a) 
6692          @staticmethod 
6693 -        def vequals_approximately(a, b, eps=1e-12): 
6694              a, b = array(a), array(b) 
6695              idx = np.where(a < 0)[0] 
6696              if len(idx): 
6697                  a[idx], b[idx] = -1 * a[idx], -1 * b[idx] 
6698              return (np.all(a - eps < b) and np.all(b < a + eps) 
6699                      ) or (np.all((1 - eps) * a < b) and np.all(b < (1 + eps) * a)) 
6700          @staticmethod 
6701 -        def expms(A, eig=np.linalg.eigh): 
6702              """matrix exponential for a symmetric matrix""" 
6703              # TODO: check that this works reliably for low rank matrices 
6704              # first: symmetrize A 
6705              D, B = eig(A) 
6706              return np.dot(B, (np.exp(D) * B).T) 
6707          @staticmethod 
6708 -        def amax(vec, vec_or_scalar): 
6709              return array(Misc.MathHelperFunctions.max(vec, vec_or_scalar)) 
6710          @staticmethod 
6711 -        def max(vec, vec_or_scalar): 
6712              b = vec_or_scalar 
6713              if np.isscalar(b): 
6714                  m = [max(x, b) for x in vec] 
6715              else: 
6716                  m = [max(vec[i], b[i]) for i in rglen((vec))] 
6717              return m 
6718          @staticmethod 
6719 -        def minmax(val, min_val, max_val): 
6720              assert min_val <= max_val 
6721              return min((max_val, max((val, min_val)))) 
6722          @staticmethod 
6723 -        def amin(vec_or_scalar, vec_or_scalar2): 
6724              return array(Misc.MathHelperFunctions.min(vec_or_scalar, vec_or_scalar2)) 
6725          @staticmethod 
6726 -        def min(a, b): 
6727              iss = np.isscalar 
6728              if iss(a) and iss(b): 
6729                  return min(a, b) 
6730              if iss(a): 
6731                  a, b = b, a 
6732              # now only b can be still a scalar 
6733              if iss(b): 
6734                  return [min(x, b) for x in a] 
6735              else:  # two non-scalars must have the same length 
6736                  return [min(a[i], b[i]) for i in rglen((a))] 
6737          @staticmethod 
6738 -        def norm(vec, expo=2): 
6739              return sum(vec**expo)**(1 / expo) 
6740          @staticmethod 
6741 -        def apos(x, lower=0): 
6742              """clips argument (scalar or array) from below at lower""" 
6743              if lower == 0: 
6744                  return (x > 0) * x 
6745              else: 
6746                  return lower + (x > lower) * (x - lower) 
6747          @staticmethod 
6748 -        def prctile(data, p_vals=[0, 25, 50, 75, 100], sorted_=False): 
6749              """``prctile(data, 50)`` returns the median, but p_vals can 
6750              also be a sequence. 
6751   
6752              Provides for small samples better values than matplotlib.mlab.prctile, 
6753              however also slower. 
6754   
6755              """ 
6756              ps = [p_vals] if np.isscalar(p_vals) else p_vals 
6757   
6758              if not sorted_: 
6759                  data = sorted(data) 
6760              n = len(data) 
6761              d = [] 
6762              for p in ps: 
6763                  fi = p * n / 100 - 0.5 
6764                  if fi <= 0:  # maybe extrapolate? 
6765                      d.append(data[0]) 
6766                  elif fi >= n - 1: 
6767                      d.append(data[-1]) 
6768                  else: 
6769                      i = int(fi) 
6770                      d.append((i + 1 - fi) * data[i] + (fi - i) * data[i + 1]) 
6771              return d[0] if np.isscalar(p_vals) else d 
6772          @staticmethod 
6773 -        def sround(nb):  # TODO: to be vectorized 
6774              """return stochastic round: floor(nb) + (rand()<remainder(nb))""" 
6775              return nb // 1 + (np.random.rand(1)[0] < (nb % 1)) 
6776   
6777          @staticmethod 
6778 -        def cauchy_with_variance_one(): 
6779              n = np.random.randn() / np.random.randn() 
6780              while abs(n) > 1000: 
6781                  n = np.random.randn() / np.random.randn() 
6782              return n / 25 
6783          @staticmethod 
6784 -        def standard_finite_cauchy(size=1): 
6785              try: 
6786                  l = len(size) 
6787              except TypeError: 
6788                  l = 0 
6789   
6790              if l == 0: 
6791                  return array([Mh.cauchy_with_variance_one() for _i in xrange(size)]) 
6792              elif l == 1: 
6793                  return array([Mh.cauchy_with_variance_one() for _i in xrange(size[0])]) 
6794              elif l == 2: 
6795                  return array([[Mh.cauchy_with_variance_one() for _i in xrange(size[1])] 
6796                               for _j in xrange(size[0])]) 
6797              else: 
6798                  raise _Error('len(size) cannot be large than two') 
6799   
6800   
6801      @staticmethod 
6802 -    def likelihood(x, m=None, Cinv=None, sigma=1, detC=None): 
6803          """return likelihood of x for the normal density N(m, sigma**2 * Cinv**-1)""" 
6804          # testing: MC integrate must be one: mean(p(x_i)) * volume(where x_i are uniformely sampled) 
6805          # for i in range(3): print mean([cma.likelihood(20*r-10, dim * [0], None, 3) for r in rand(10000,dim)]) * 20**dim 
6806          if m is None: 
6807              dx = x 
6808          else: 
6809              dx = x - m  # array(x) - array(m) 
6810          n = len(x) 
6811          s2pi = (2 * np.pi)**(n / 2.) 
6812          if Cinv is None: 
6813              return exp(-sum(dx**2) / sigma**2 / 2) / s2pi / sigma**n 
6814          if detC is None: 
6815              detC = 1. / np.linalg.linalg.det(Cinv) 
6816          return  exp(-np.dot(dx, np.dot(Cinv, dx)) / sigma**2 / 2) / s2pi / abs(detC)**0.5 / sigma**n 
6817   
6818      @staticmethod 
6819 -    def loglikelihood(self, x, previous=False): 
6820          """return log-likelihood of `x` regarding the current sample distribution""" 
6821          # testing of original fct: MC integrate must be one: mean(p(x_i)) * volume(where x_i are uniformely sampled) 
6822          # for i in range(3): print mean([cma.likelihood(20*r-10, dim * [0], None, 3) for r in rand(10000,dim)]) * 20**dim 
6823          # TODO: test this!! 
6824          # c=cma.fmin... 
6825          # c[3]['cma'].loglikelihood(...) 
6826   
6827          if previous and hasattr(self, 'lastiter'): 
6828              sigma = self.lastiter.sigma 
6829              Crootinv = self.lastiter._Crootinv 
6830              xmean = self.lastiter.mean 
6831              D = self.lastiter.D 
6832          elif previous and self.countiter > 1: 
6833              raise _Error('no previous distribution parameters stored, check options importance_mixing') 
6834          else: 
6835              sigma = self.sigma 
6836              Crootinv = self._Crootinv 
6837              xmean = self.mean 
6838              D = self.D 
6839   
6840          dx = array(x) - xmean  # array(x) - array(m) 
6841          n = self.N 
6842          logs2pi = n * log(2 * np.pi) / 2. 
6843          logdetC = 2 * sum(log(D)) 
6844          dx = np.dot(Crootinv, dx) 
6845          res = -sum(dx**2) / sigma**2 / 2 - logs2pi - logdetC / 2 - n * log(sigma) 
6846          if 1 < 3:  # testing 
6847              s2pi = (2 * np.pi)**(n / 2.) 
6848              detC = np.prod(D)**2 
6849              res2 = -sum(dx**2) / sigma**2 / 2 - log(s2pi * abs(detC)**0.5 * sigma**n) 
6850              assert res2 < res + 1e-8 or res2 > res - 1e-8 
6851          return res 
6852   
6853      # ____________________________________________________________ 
6854      # ____________________________________________________________ 
6855      # 
6856      # C and B are arrays rather than matrices, because they are 
6857      # addressed via B[i][j], matrices can only be addressed via B[i,j] 
6858   
6859      # tred2(N, B, diagD, offdiag); 
6860      # tql2(N, diagD, offdiag, B); 
6861   
6862   
6863      # Symmetric Householder reduction to tridiagonal form, translated from JAMA package. 
6864      @staticmethod 
6865 -    def eig(C): 
6866          """eigendecomposition of a symmetric matrix, much slower than 
6867          `numpy.linalg.eigh`, return ``(EVals, Basis)``, the eigenvalues 
6868          and an orthonormal basis of the corresponding eigenvectors, where 
6869   
6870              ``Basis[i]`` 
6871                  the i-th row of ``Basis`` 
6872              columns of ``Basis``, ``[Basis[j][i] for j in range(len(Basis))]`` 
6873                  the i-th eigenvector with eigenvalue ``EVals[i]`` 
6874   
6875          """ 
6876   
6877      # class eig(object): 
6878      #     def __call__(self, C): 
6879   
6880      # Householder transformation of a symmetric matrix V into tridiagonal form. 
6881          # -> n             : dimension 
6882          # -> V             : symmetric nxn-matrix 
6883          # <- V             : orthogonal transformation matrix: 
6884          #                    tridiag matrix == V * V_in * V^t 
6885          # <- d             : diagonal 
6886          # <- e[0..n-1]     : off diagonal (elements 1..n-1) 
6887   
6888          # Symmetric tridiagonal QL algorithm, iterative 
6889          # Computes the eigensystem from a tridiagonal matrix in roughtly 3N^3 operations 
6890          # -> n     : Dimension. 
6891          # -> d     : Diagonale of tridiagonal matrix. 
6892          # -> e[1..n-1] : off-diagonal, output from Householder 
6893          # -> V     : matrix output von Householder 
6894          # <- d     : eigenvalues 
6895          # <- e     : garbage? 
6896          # <- V     : basis of eigenvectors, according to d 
6897   
6898   
6899          #  tred2(N, B, diagD, offdiag); B=C on input 
6900          #  tql2(N, diagD, offdiag, B); 
6901   
6902          #  private void tred2 (int n, double V[][], double d[], double e[]) { 
6903          def tred2 (n, V, d, e): 
6904              #  This is derived from the Algol procedures tred2 by 
6905              #  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for 
6906              #  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding 
6907              #  Fortran subroutine in EISPACK. 
6908   
6909              num_opt = False  # factor 1.5 in 30-D 
6910   
6911              for j in range(n): 
6912                  d[j] = V[n - 1][j]  # d is output argument 
6913   
6914              # Householder reduction to tridiagonal form. 
6915   
6916              for i in range(n - 1, 0, -1): 
6917                  # Scale to avoid under/overflow. 
6918                  h = 0.0 
6919                  if not num_opt: 
6920                      scale = 0.0 
6921                      for k in range(i): 
6922                          scale = scale + abs(d[k]) 
6923                  else: 
6924                      scale = sum(abs(d[0:i])) 
6925   
6926                  if scale == 0.0: 
6927                      e[i] = d[i - 1] 
6928                      for j in range(i): 
6929                          d[j] = V[i - 1][j] 
6930                          V[i][j] = 0.0 
6931                          V[j][i] = 0.0 
6932                  else: 
6933   
6934                      # Generate Householder vector. 
6935                      if not num_opt: 
6936                          for k in range(i): 
6937                              d[k] /= scale 
6938                              h += d[k] * d[k] 
6939                      else: 
6940                          d[:i] /= scale 
6941                          h = np.dot(d[:i], d[:i]) 
6942   
6943                      f = d[i - 1] 
6944                      g = h**0.5 
6945   
6946                      if f > 0: 
6947                          g = -g 
6948   
6949                      e[i] = scale * g 
6950                      h = h - f * g 
6951                      d[i - 1] = f - g 
6952                      if not num_opt: 
6953                          for j in range(i): 
6954                              e[j] = 0.0 
6955                      else: 
6956                          e[:i] = 0.0 
6957   
6958                      # Apply similarity transformation to remaining columns. 
6959   
6960                      for j in range(i): 
6961                          f = d[j] 
6962                          V[j][i] = f 
6963                          g = e[j] + V[j][j] * f 
6964                          if not num_opt: 
6965                              for k in range(j + 1, i): 
6966                                  g += V[k][j] * d[k] 
6967                                  e[k] += V[k][j] * f 
6968                              e[j] = g 
6969                          else: 
6970                              e[j + 1:i] += V.T[j][j + 1:i] * f 
6971                              e[j] = g + np.dot(V.T[j][j + 1:i], d[j + 1:i]) 
6972   
6973                      f = 0.0 
6974                      if not num_opt: 
6975                          for j in range(i): 
6976                              e[j] /= h 
6977                              f += e[j] * d[j] 
6978                      else: 
6979                          e[:i] /= h 
6980                          f += np.dot(e[:i], d[:i]) 
6981   
6982                      hh = f / (h + h) 
6983                      if not num_opt: 
6984                          for j in range(i): 
6985                              e[j] -= hh * d[j] 
6986                      else: 
6987                          e[:i] -= hh * d[:i] 
6988   
6989                      for j in range(i): 
6990                          f = d[j] 
6991                          g = e[j] 
6992                          if not num_opt: 
6993                              for k in range(j, i): 
6994                                  V[k][j] -= (f * e[k] + g * d[k]) 
6995                          else: 
6996                              V.T[j][j:i] -= (f * e[j:i] + g * d[j:i]) 
6997   
6998                          d[j] = V[i - 1][j] 
6999                          V[i][j] = 0.0 
7000   
7001                  d[i] = h 
7002              # end for i-- 
7003   
7004              # Accumulate transformations. 
7005   
7006              for i in range(n - 1): 
7007                  V[n - 1][i] = V[i][i] 
7008                  V[i][i] = 1.0 
7009                  h = d[i + 1] 
7010                  if h != 0.0: 
7011                      if not num_opt: 
7012                          for k in range(i + 1): 
7013                              d[k] = V[k][i + 1] / h 
7014                      else: 
7015                          d[:i + 1] = V.T[i + 1][:i + 1] / h 
7016   
7017                      for j in range(i + 1): 
7018                          if not num_opt: 
7019                              g = 0.0 
7020                              for k in range(i + 1): 
7021                                  g += V[k][i + 1] * V[k][j] 
7022                              for k in range(i + 1): 
7023                                  V[k][j] -= g * d[k] 
7024                          else: 
7025                              g = np.dot(V.T[i + 1][0:i + 1], V.T[j][0:i + 1]) 
7026                              V.T[j][:i + 1] -= g * d[:i + 1] 
7027   
7028                  if not num_opt: 
7029                      for k in range(i + 1): 
7030                          V[k][i + 1] = 0.0 
7031                  else: 
7032                      V.T[i + 1][:i + 1] = 0.0 
7033   
7034   
7035              if not num_opt: 
7036                  for j in range(n): 
7037                      d[j] = V[n - 1][j] 
7038                      V[n - 1][j] = 0.0 
7039              else: 
7040                  d[:n] = V[n - 1][:n] 
7041                  V[n - 1][:n] = 0.0 
7042   
7043              V[n - 1][n - 1] = 1.0 
7044              e[0] = 0.0 
7045   
7046   
7047          # Symmetric tridiagonal QL algorithm, taken from JAMA package. 
7048          # private void tql2 (int n, double d[], double e[], double V[][]) { 
7049          # needs roughly 3N^3 operations 
7050          def tql2 (n, d, e, V): 
7051   
7052              #  This is derived from the Algol procedures tql2, by 
7053              #  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for 
7054              #  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding 
7055              #  Fortran subroutine in EISPACK. 
7056   
7057              num_opt = False  # using vectors from numpy makes it faster 
7058   
7059              if not num_opt: 
7060                  for i in range(1, n):  # (int i = 1; i < n; i++): 
7061                      e[i - 1] = e[i] 
7062              else: 
7063                  e[0:n - 1] = e[1:n] 
7064              e[n - 1] = 0.0 
7065   
7066              f = 0.0 
7067              tst1 = 0.0 
7068              eps = 2.0**-52.0 
7069              for l in range(n):  # (int l = 0; l < n; l++) { 
7070   
7071                  # Find small subdiagonal element 
7072   
7073                  tst1 = max(tst1, abs(d[l]) + abs(e[l])) 
7074                  m = l 
7075                  while m < n: 
7076                      if abs(e[m]) <= eps * tst1: 
7077                          break 
7078                      m += 1 
7079   
7080                  # If m == l, d[l] is an eigenvalue, 
7081                  # otherwise, iterate. 
7082   
7083                  if m > l: 
7084                      iiter = 0 
7085                      while 1:  # do { 
7086                          iiter += 1  # (Could check iteration count here.) 
7087   
7088                          # Compute implicit shift 
7089   
7090                          g = d[l] 
7091                          p = (d[l + 1] - g) / (2.0 * e[l]) 
7092                          r = (p**2 + 1)**0.5  # hypot(p,1.0) 
7093                          if p < 0: 
7094                              r = -r 
7095   
7096                          d[l] = e[l] / (p + r) 
7097                          d[l + 1] = e[l] * (p + r) 
7098                          dl1 = d[l + 1] 
7099                          h = g - d[l] 
7100                          if not num_opt: 
7101                              for i in range(l + 2, n): 
7102                                  d[i] -= h 
7103                          else: 
7104                              d[l + 2:n] -= h 
7105   
7106                          f = f + h 
7107   
7108                          # Implicit QL transformation. 
7109   
7110                          p = d[m] 
7111                          c = 1.0 
7112                          c2 = c 
7113                          c3 = c 
7114                          el1 = e[l + 1] 
7115                          s = 0.0 
7116                          s2 = 0.0 
7117   
7118                          # hh = V.T[0].copy()  # only with num_opt 
7119                          for i in range(m - 1, l - 1, -1):  # (int i = m-1; i >= l; i--) { 
7120                              c3 = c2 
7121                              c2 = c 
7122                              s2 = s 
7123                              g = c * e[i] 
7124                              h = c * p 
7125                              r = (p**2 + e[i]**2)**0.5  # hypot(p,e[i]) 
7126                              e[i + 1] = s * r 
7127                              s = e[i] / r 
7128                              c = p / r 
7129                              p = c * d[i] - s * g 
7130                              d[i + 1] = h + s * (c * g + s * d[i]) 
7131   
7132                              # Accumulate transformation. 
7133   
7134                              if not num_opt:  # overall factor 3 in 30-D 
7135                                  for k in range(n):  # (int k = 0; k < n; k++) { 
7136                                      h = V[k][i + 1] 
7137                                      V[k][i + 1] = s * V[k][i] + c * h 
7138                                      V[k][i] = c * V[k][i] - s * h 
7139                              else:  # about 20% faster in 10-D 
7140                                  hh = V.T[i + 1].copy() 
7141                                  # hh[:] = V.T[i+1][:] 
7142                                  V.T[i + 1] = s * V.T[i] + c * hh 
7143                                  V.T[i] = c * V.T[i] - s * hh 
7144                                  # V.T[i] *= c 
7145                                  # V.T[i] -= s * hh 
7146   
7147                          p = -s * s2 * c3 * el1 * e[l] / dl1 
7148                          e[l] = s * p 
7149                          d[l] = c * p 
7150   
7151                          # Check for convergence. 
7152                          if abs(e[l]) <= eps * tst1: 
7153                              break 
7154                      # } while (Math.abs(e[l]) > eps*tst1); 
7155   
7156                  d[l] = d[l] + f 
7157                  e[l] = 0.0 
7158   
7159   
7160              # Sort eigenvalues and corresponding vectors. 
7161          # tql2 
7162   
7163          N = len(C[0]) 
7164          if 1 < 3: 
7165              V = [[x[i] for i in xrange(N)] for x in C]  # copy each "row" 
7166              d = N * [0.] 
7167              e = N * [0.] 
7168   
7169          tred2(N, V, d, e) 
7170          tql2(N, d, e, V) 
7171          return (array(d), array(V)) 
7172  Mh = Misc.MathHelperFunctions 
7173   
7174  # from new_stuff import * 
7175   
7176 -def pprint(to_be_printed): 
7177      """nicely formated print""" 
7178      try: 
7179          import pprint as pp 
7180          # generate an instance PrettyPrinter 
7181          # pp.PrettyPrinter().pprint(to_be_printed) 
7182          pp.pprint(to_be_printed) 
7183      except ImportError: 
7184          if isinstance(to_be_printed, dict): 
7185              print('{') 
7186              for k, v in to_be_printed.items(): 
7187                  print("'" + k + "'" if isinstance(k, basestring) else k, 
7188                        ': ', 
7189                        "'" + v + "'" if isinstance(k, basestring) else v, 
7190                        sep="") 
7191              print('}') 
7192          else: 
7193              print('could not import pprint module, will apply regular print') 
7194              print(to_be_printed) 
7195   
7196  pp = pprint 
7197   
7198 -class Rotation(object): 
7199      """Rotation class that implements an orthogonal linear transformation, 
7200      one for each dimension. Used to implement non-separable test functions. 
7201   
7202      Example: 
7203   
7204      >>> import cma, numpy as np 
7205      >>> R = cma.Rotation() 
7206      >>> R2 = cma.Rotation() # another rotation 
7207      >>> x = np.array((1,2,3)) 
7208      >>> print(R(R(x), inverse=1)) 
7209      [ 1.  2.  3.] 
7210   
7211      """ 
7212      dicMatrices = {}  # store matrix if necessary, for each dimension 
7213 -    def __init__(self): 
7214          self.dicMatrices = {}  # otherwise there might be shared bases which is probably not what we want 
7215 -    def __call__(self, x, inverse=False):  # function when calling an object 
7216          """Rotates the input array `x` with a fixed rotation matrix 
7217             (``self.dicMatrices['str(len(x))']``) 
7218          """ 
7219          N = x.shape[0]  # can be an array or matrix, TODO: accept also a list of arrays? 
7220          if str(N) not in self.dicMatrices:  # create new N-basis for once and all 
7221              B = np.random.randn(N, N) 
7222              for i in xrange(N): 
7223                  for j in xrange(0, i): 
7224                      B[i] -= np.dot(B[i], B[j]) * B[j] 
7225                  B[i] /= sum(B[i]**2)**0.5 
7226              self.dicMatrices[str(N)] = B 
7227          if inverse: 
7228              return np.dot(self.dicMatrices[str(N)].T, x)  # compute rotation 
7229          else: 
7230              return np.dot(self.dicMatrices[str(N)], x)  # compute rotation 
7231  # Use rotate(x) to rotate x 
7232  rotate = Rotation() 
7233   
7234  # ____________________________________________________________ 
7235  # ____________________________________________________________ 
7236  # 
7237 -class FitnessFunctions(object): 
7238      """ versatile container for test objective functions """ 
7239   
7240 -    def __init__(self): 
7241          self.counter = 0  # number of calls or any other practical use 
7242 -    def rot(self, x, fun, rot=1, args=()): 
7243          """returns ``fun(rotation(x), *args)``, ie. `fun` applied to a rotated argument""" 
7244          if len(np.shape(array(x))) > 1:  # parallelized 
7245              res = [] 
7246              for x in x: 
7247                  res.append(self.rot(x, fun, rot, args)) 
7248              return res 
7249   
7250          if rot: 
7251              return fun(rotate(x, *args)) 
7252          else: 
7253              return fun(x) 
7254 -    def somenan(self, x, fun, p=0.1): 
7255          """returns sometimes np.NaN, otherwise fun(x)""" 
7256          if np.random.rand(1) < p: 
7257              return np.NaN 
7258          else: 
7259              return fun(x) 
7260 -    def rand(self, x): 
7261          """Random test objective function""" 
7262          return np.random.random(1)[0] 
7263 -    def linear(self, x): 
7264          return -x[0] 
7265 -    def lineard(self, x): 
7266          if 1 < 3 and any(array(x) < 0): 
7267              return np.nan 
7268          if 1 < 3 and sum([ (10 + i) * x[i] for i in rglen(x)]) > 50e3: 
7269              return np.nan 
7270          return -sum(x) 
7271 -    def sphere(self, x): 
7272          """Sphere (squared norm) test objective function""" 
7273          # return np.random.rand(1)[0]**0 * sum(x**2) + 1 * np.random.rand(1)[0] 
7274          return sum((x + 0)**2) 
7275 -    def sphere_pos(self, x): 
7276          """Sphere (squared norm) test objective function""" 
7277          # return np.random.rand(1)[0]**0 * sum(x**2) + 1 * np.random.rand(1)[0] 
7278          c = 0.0 
7279          if x[0] < c: 
7280              return np.nan 
7281          return -c**2 + sum((x + 0)**2) 
7282 -    def spherewithoneconstraint(self, x): 
7283          return sum((x + 0)**2) if x[0] > 1 else np.nan 
7284 -    def elliwithoneconstraint(self, x, idx=[-1]): 
7285          return self.ellirot(x) if all(array(x)[idx] > 1) else np.nan 
7286   
7287 -    def spherewithnconstraints(self, x): 
7288          return sum((x + 0)**2) if all(array(x) > 1) else np.nan 
7289      # zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz 
7290 -    def noisysphere(self, x, noise=2.10e-9, cond=1.0, noise_offset=0.10): 
7291          """noise=10 does not work with default popsize, noise handling does not help """ 
7292          return self.elli(x, cond=cond) * (1 + noise * np.random.randn() / len(x)) + noise_offset * np.random.rand() 
7293 -    def spherew(self, x): 
7294          """Sphere (squared norm) with sum x_i = 1 test objective function""" 
7295          # return np.random.rand(1)[0]**0 * sum(x**2) + 1 * np.random.rand(1)[0] 
7296          # s = sum(abs(x)) 
7297          # return sum((x/s+0)**2) - 1/len(x) 
7298          # return sum((x/s)**2) - 1/len(x) 
7299          return -0.01 * x[0] + abs(x[0])**-2 * sum(x[1:]**2) 
7300 -    def partsphere(self, x): 
7301          """Sphere (squared norm) test objective function""" 
7302          self.counter += 1 
7303          # return np.random.rand(1)[0]**0 * sum(x**2) + 1 * np.random.rand(1)[0] 
7304          dim = len(x) 
7305          x = array([x[i % dim] for i in range(2 * dim)]) 
7306          N = 8 
7307          i = self.counter % dim 
7308          # f = sum(x[i:i + N]**2) 
7309          f = sum(x[np.random.randint(dim, size=N)]**2) 
7310          return f 
7311 -    def sectorsphere(self, x): 
7312          """asymmetric Sphere (squared norm) test objective function""" 
7313          return sum(x**2) + (1e6 - 1) * sum(x[x < 0]**2) 
7314 -    def cornersphere(self, x): 
7315          """Sphere (squared norm) test objective function constraint to the corner""" 
7316          nconstr = len(x) - 0 
7317          if any(x[:nconstr] < 1): 
7318              return np.NaN 
7319          return sum(x**2) - nconstr 
7320 -    def cornerelli(self, x): 
7321          """ """ 
7322          if any(x < 1): 
7323              return np.NaN 
7324          return self.elli(x) - self.elli(np.ones(len(x))) 
7325 -    def cornerellirot(self, x): 
7326          """ """ 
7327          if any(x < 1): 
7328              return np.NaN 
7329          return self.ellirot(x) 
7330 -    def normalSkew(self, f): 
7331          N = np.random.randn(1)[0]**2 
7332          if N < 1: 
7333              N = f * N  # diminish blow up lower part 
7334          return N 
7335 -    def noiseC(self, x, func=sphere, fac=10, expon=0.8): 
7336          f = func(self, x) 
7337          N = np.random.randn(1)[0] / np.random.randn(1)[0] 
7338          return max(1e-19, f + (float(fac) / len(x)) * f**expon * N) 
7339 -    def noise(self, x, func=sphere, fac=10, expon=1): 
7340          f = func(self, x) 
7341          # R = np.random.randn(1)[0] 
7342          R = np.log10(f) + expon * abs(10 - np.log10(f)) * np.random.rand(1)[0] 
7343          # sig = float(fac)/float(len(x)) 
7344          # R = log(f) + 0.5*log(f) * random.randn(1)[0] 
7345          # return max(1e-19, f + sig * (f**np.log10(f)) * np.exp(R)) 
7346          # return max(1e-19, f * np.exp(sig * N / f**expon)) 
7347          # return max(1e-19, f * normalSkew(f**expon)**sig) 
7348          return f + 10**R  # == f + f**(1+0.5*RN) 
7349 -    def cigar(self, x, rot=0, cond=1e6, noise=0): 
7350          """Cigar test objective function""" 
7351          if rot: 
7352              x = rotate(x) 
7353          x = [x] if np.isscalar(x[0]) else x  # scalar into list 
7354          f = [(x[0]**2 + cond * sum(x[1:]**2)) * np.exp(noise * np.random.randn(1)[0] / len(x)) for x in x] 
7355          return f if len(f) > 1 else f[0]  # 1-element-list into scalar 
7356 -    def tablet(self, x, rot=0): 
7357          """Tablet test objective function""" 
7358          if rot: 
7359              x = rotate(x) 
7360          x = [x] if np.isscalar(x[0]) else x  # scalar into list 
7361          f = [1e6 * x[0]**2 + sum(x[1:]**2) for x in x] 
7362          return f if len(f) > 1 else f[0]  # 1-element-list into scalar 
7363 -    def cigtab(self, y): 
7364          """Cigtab test objective function""" 
7365          X = [y] if np.isscalar(y[0]) else y 
7366          f = [1e-4 * x[0]**2 + 1e4 * x[1]**2 + sum(x[2:]**2) for x in X] 
7367          return f if len(f) > 1 else f[0] 
7368 -    def twoaxes(self, y): 
7369          """Cigtab test objective function""" 
7370          X = [y] if np.isscalar(y[0]) else y 
7371          N2 = len(X[0]) // 2 
7372          f = [1e6 * sum(x[0:N2]**2) + sum(x[N2:]**2) for x in X] 
7373          return f if len(f) > 1 else f[0] 
7374 -    def ellirot(self, x): 
7375          return fcts.elli(array(x), 1) 
7376 -    def hyperelli(self, x): 
7377          N = len(x) 
7378          return sum((np.arange(1, N + 1) * x)**2) 
7379 -    def elli(self, x, rot=0, xoffset=0, cond=1e6, actuator_noise=0.0, both=False): 
7380          """Ellipsoid test objective function""" 
7381          if not np.isscalar(x[0]):  # parallel evaluation 
7382              return [self.elli(xi, rot) for xi in x]  # could save 20% overall 
7383          if rot: 
7384              x = rotate(x) 
7385          N = len(x) 
7386          if actuator_noise: 
7387              x = x + actuator_noise * np.random.randn(N) 
7388   
7389          ftrue = sum(cond**(np.arange(N) / (N - 1.)) * (x + xoffset)**2) 
7390   
7391          alpha = 0.49 + 1. / N 
7392          beta = 1 
7393          felli = np.random.rand(1)[0]**beta * ftrue * \ 
7394                  max(1, (10.**9 / (ftrue + 1e-99))**(alpha * np.random.rand(1)[0])) 
7395          # felli = ftrue + 1*np.random.randn(1)[0] / (1e-30 + 
7396          #                                           np.abs(np.random.randn(1)[0]))**0 
7397          if both: 
7398              return (felli, ftrue) 
7399          else: 
7400              # return felli  # possibly noisy value 
7401              return ftrue  # + np.random.randn() 
7402 -    def elliconstraint(self, x, cfac=1e8, tough=True, cond=1e6): 
7403          """ellipsoid test objective function with "constraints" """ 
7404          N = len(x) 
7405          f = sum(cond**(np.arange(N)[-1::-1] / (N - 1)) * x**2) 
7406          cvals = (x[0] + 1, 
7407                   x[0] + 1 + 100 * x[1], 
7408                   x[0] + 1 - 100 * x[1]) 
7409          if tough: 
7410              f += cfac * sum(max(0, c) for c in cvals) 
7411          else: 
7412              f += cfac * sum(max(0, c + 1e-3)**2 for c in cvals) 
7413          return f 
7414 -    def rosen(self, x, alpha=1e2): 
7415          """Rosenbrock test objective function""" 
7416          x = [x] if np.isscalar(x[0]) else x  # scalar into list 
7417          f = [sum(alpha * (x[:-1]**2 - x[1:])**2 + (1. - x[:-1])**2) for x in x] 
7418          return f if len(f) > 1 else f[0]  # 1-element-list into scalar 
7419 -    def diffpow(self, x, rot=0): 
7420          """Diffpow test objective function""" 
7421          N = len(x) 
7422          if rot: 
7423              x = rotate(x) 
7424          return sum(np.abs(x)**(2. + 4.*np.arange(N) / (N - 1.)))**0.5 
7425 -    def rosenelli(self, x): 
7426          N = len(x) 
7427          return self.rosen(x[:N / 2]) + self.elli(x[N / 2:], cond=1) 
7428 -    def ridge(self, x, expo=2): 
7429          x = [x] if np.isscalar(x[0]) else x  # scalar into list 
7430          f = [x[0] + 100 * np.sum(x[1:]**2)**(expo / 2.) for x in x] 
7431          return f if len(f) > 1 else f[0]  # 1-element-list into scalar 
7432 -    def ridgecircle(self, x, expo=0.5): 
7433          """happy cat by HG Beyer""" 
7434          a = len(x) 
7435          s = sum(x**2) 
7436          return ((s - a)**2)**(expo / 2) + s / a + sum(x) / a 
7437 -    def happycat(self, x, alpha=1. / 8): 
7438          s = sum(x**2) 
7439          return ((s - len(x))**2)**alpha + (s / 2 + sum(x)) / len(x) + 0.5 
7440 -    def flat(self, x): 
7441          return 1 
7442          return 1 if np.random.rand(1) < 0.9 else 1.1 
7443          return np.random.randint(1, 30) 
7444 -    def branin(self, x): 
7445          # in [0,15]**2 
7446          y = x[1] 
7447          x = x[0] + 5 
7448          return (y - 5.1 * x**2 / 4 / np.pi**2 + 5 * x / np.pi - 6)**2 + 10 * (1 - 1 / 8 / np.pi) * np.cos(x) + 10 - 0.397887357729738160000 
7449 -    def goldsteinprice(self, x): 
7450          x1 = x[0] 
7451          x2 = x[1] 
7452          return (1 + (x1 + x2 + 1)**2 * (19 - 14 * x1 + 3 * x1**2 - 14 * x2 + 6 * x1 * x2 + 3 * x2**2)) * ( 
7453                  30 + (2 * x1 - 3 * x2)**2 * (18 - 32 * x1 + 12 * x1**2 + 48 * x2 - 36 * x1 * x2 + 27 * x2**2)) - 3 
7454 -    def griewank(self, x): 
7455          # was in [-600 600] 
7456          x = (600. / 5) * x 
7457          return 1 - np.prod(np.cos(x / sqrt(1. + np.arange(len(x))))) + sum(x**2) / 4e3 
7458 -    def rastrigin(self, x): 
7459          """Rastrigin test objective function""" 
7460          if not np.isscalar(x[0]): 
7461              N = len(x[0]) 
7462              return [10 * N + sum(xi**2 - 10 * np.cos(2 * np.pi * xi)) for xi in x] 
7463              # return 10*N + sum(x**2 - 10*np.cos(2*np.pi*x), axis=1) 
7464          N = len(x) 
7465          return 10 * N + sum(x**2 - 10 * np.cos(2 * np.pi * x)) 
7466 -    def schaffer(self, x): 
7467          """ Schaffer function x0 in [-100..100]""" 
7468          N = len(x) 
7469          s = x[0:N - 1]**2 + x[1:N]**2 
7470          return sum(s**0.25 * (np.sin(50 * s**0.1)**2 + 1)) 
7471   
7472 -    def schwefelelli(self, x): 
7473          s = 0 
7474          f = 0 
7475          for i in rglen(x): 
7476              s += x[i] 
7477              f += s**2 
7478          return f 
7479 -    def schwefelmult(self, x, pen_fac=1e4): 
7480          """multimodal Schwefel function with domain -500..500""" 
7481          y = [x] if np.isscalar(x[0]) else x 
7482          N = len(y[0]) 
7483          f = array([418.9829 * N - 1.27275661e-5 * N - sum(x * np.sin(np.abs(x)**0.5)) 
7484                  + pen_fac * sum((abs(x) > 500) * (abs(x) - 500)**2) for x in y]) 
7485          return f if len(f) > 1 else f[0] 
7486 -    def optprob(self, x): 
7487          n = np.arange(len(x)) + 1 
7488          f = n * x * (1 - x)**(n - 1) 
7489          return sum(1 - f) 
7490 -    def lincon(self, x, theta=0.01): 
7491          """ridge like linear function with one linear constraint""" 
7492          if x[0] < 0: 
7493              return np.NaN 
7494          return theta * x[1] + x[0] 
7495 -    def rosen_nesterov(self, x, rho=100): 
7496          """needs exponential number of steps in a non-increasing f-sequence. 
7497   
7498          x_0 = (-1,1,...,1) 
7499          See Jarre (2011) "On Nesterov's Smooth Chebyshev-Rosenbrock Function" 
7500   
7501          """ 
7502          f = 0.25 * (x[0] - 1)**2 
7503          f += rho * sum((x[1:] - 2 * x[:-1]**2 + 1)**2) 
7504          return f 
7505 -    def powel_singular(self, x): 
7506          # ((8 * np.sin(7 * (x[i] - 0.9)**2)**2 ) + (6 * np.sin())) 
7507          res = np.sum((x[i - 1] + 10 * x[i])**2 + 5 * (x[i + 1] - x[i + 2])**2 + 
7508                       (x[i] - 2 * x[i + 1])**4 + 10 * (x[i - 1] - x[i + 2])**4 
7509                       for i in xrange(1, len(x) - 2)) 
7510          return 1 + res 
7511 -    def styblinski_tang(self, x): 
7512          return (39.1661657037714171054273576010019 * len(x))**1 + sum(x**4 - 16*x**2 + 5*x) / 2 
7513 -    def bukin(self, x): 
7514          """Bukin function from Wikipedia, generalized simplistically from 2-D. 
7515   
7516          http://en.wikipedia.org/wiki/Test_functions_for_optimization""" 
7517          s = 0 
7518          for k in range((1+len(x)) // 2): 
7519              z = x[2 * k] 
7520              y = x[min((2*k + 1, len(x)-1))] 
7521              s += 100 * np.abs(y - 0.01 * z**2)**0.5 + 0.01 * np.abs(z + 10) 
7522          return s 
7523   
7524  fcts = FitnessFunctions() 
7525  Fcts = fcts  # for cross compatibility, as if the functions were static members of class Fcts 
7526 -def felli(x): 
7527      """unbound test function, needed to test multiprocessor""" 
7528      return sum(1e6**(np.arange(len(x)) / (len(x) - 1)) * (np.array(x, copy=False))**2) 
7529   
7530   
7531  # ____________________________________________ 
7532  # ____________________________________________________________ 
7533 -def _test(module=None):  # None is fine when called from inside the module 
7534      import doctest 
7535      print(doctest.testmod(module))  # this is pretty coool! 
7536 -def process_doctest_output(stream=None): 
7537      """ """ 
7538      import fileinput 
7539      s1 = "" 
7540      s2 = "" 
7541      s3 = "" 
7542      state = 0 
7543      for line in fileinput.input(stream):  # takes argv as file or stdin 
7544          if 1 < 3: 
7545              s3 += line 
7546              if state < -1 and line.startswith('***'): 
7547                  print(s3) 
7548              if line.startswith('***'): 
7549                  s3 = "" 
7550   
7551          if state == -1:  # found a failed example line 
7552              s1 += '\n\n*** Failed Example:' + line 
7553              s2 += '\n\n\n'  # line 
7554              # state = 0  # wait for 'Expected:' line 
7555   
7556          if line.startswith('Expected:'): 
7557              state = 1 
7558              continue 
7559          elif line.startswith('Got:'): 
7560              state = 2 
7561              continue 
7562          elif line.startswith('***'):  # marks end of failed example 
7563              state = 0 
7564          elif line.startswith('Failed example:'): 
7565              state = -1 
7566          elif line.startswith('Exception raised'): 
7567              state = -2 
7568   
7569          # in effect more else: 
7570          if state == 1: 
7571              s1 += line + '' 
7572          if state == 2: 
7573              s2 += line + '' 
7574   
7575  # ____________________________________________________________ 
7576  # ____________________________________________________________ 
7577  # 
7578 -def main(argv=None): 
7579      """to install and/or test from the command line use:: 
7580   
7581          python cma.py [options | func dim sig0 [optkey optval][optkey optval]...] 
7582   
7583      with options being 
7584   
7585      ``--test`` (or ``-t``) to run the doctest, ``--test -v`` to get (much) verbosity. 
7586   
7587      ``install`` to install cma.py (uses setup from distutils.core). 
7588   
7589      ``--doc`` for more infos. 
7590   
7591      Or start Python or (even better) ``ipython -pylab`` and:: 
7592   
7593          import cma 
7594          cma.main('--test') 
7595          help(cma) 
7596          help(cma.fmin) 
7597          res = fmin(cma.fcts.rosen, 10 * [0], 1) 
7598          cma.plot() 
7599   
7600      Examples 
7601      ======== 
7602      Testing with the local python distribution from a command line 
7603      in a folder where ``cma.py`` can be found:: 
7604   
7605          python cma.py --test 
7606   
7607      And a single run on the Rosenbrock function:: 
7608   
7609          python cma.py rosen 10 1  # dimension initial_sigma 
7610          python cma.py plot 
7611   
7612      In the python shell:: 
7613   
7614          import cma 
7615          cma.main('--test') 
7616   
7617      """ 
7618      if argv is None: 
7619          argv = sys.argv  # should have better been sys.argv[1:] 
7620      else: 
7621          if isinstance(argv, list): 
7622              argv = ['python'] + argv  # see above 
7623          else: 
7624              argv = ['python'] + [argv] 
7625   
7626      # uncomment for unit test 
7627      # _test() 
7628      # handle input arguments, getopt might be helpful ;-) 
7629      if len(argv) >= 1:  # function and help 
7630          if len(argv) == 1 or argv[1].startswith('-h') or argv[1].startswith('--help'): 
7631              print(main.__doc__) 
7632              fun = None 
7633          elif argv[1].startswith('-t') or argv[1].startswith('--test'): 
7634              import doctest 
7635              if len(argv) > 2 and (argv[2].startswith('--v') or argv[2].startswith('-v')):  # verbose 
7636                  print('doctest for cma.py: due to different platforms and python versions') 
7637                  print('and in some cases due to a missing unique random seed') 
7638                  print('many examples will "fail". This is OK, if they give a similar') 
7639                  print('to the expected result and if no exception occurs. ') 
7640                  # if argv[1][2] == 'v': 
7641                  doctest.testmod(sys.modules[__name__], report=True)  # this is quite cool! 
7642              else:  # was: if len(argv) > 2 and (argv[2].startswith('--qu') or argv[2].startswith('-q')): 
7643                  print('doctest for cma.py: launching...') # not anymore: (it might be necessary to close the pop up window to finish) 
7644                  fn = '_cma_doctest_.txt' 
7645                  stdout = sys.stdout 
7646                  try: 
7647                      with open(fn, 'w') as f: 
7648                          sys.stdout = f 
7649                          clock = ElapsedTime() 
7650                          doctest.testmod(sys.modules[__name__], report=True)  # this is quite cool! 
7651                          t_elapsed = clock() 
7652                  finally: 
7653                      sys.stdout = stdout 
7654                  process_doctest_output(fn) 
7655                  # clean up 
7656                  try: 
7657                      import os 
7658                      for name in os.listdir('.'): 
7659                          if (name.startswith('bound_method_FitnessFunctions.rosen_of_cma.FitnessFunctions_object_at_') 
7660                              and name.endswith('.pkl')): 
7661                              os.remove(name) 
7662                  except: 
7663                      pass 
7664                  print('doctest for cma.py: finished (no other output should be seen after launching, more in file _cma_doctest_.txt)') 
7665                  print('  elapsed time [s]:', t_elapsed) 
7666              return 
7667          elif argv[1] == '--doc': 
7668              print(__doc__) 
7669              print(CMAEvolutionStrategy.__doc__) 
7670              print(fmin.__doc__) 
7671              fun = None 
7672          elif argv[1] == '--fcts': 
7673              print('List of valid function names:') 
7674              print([d for d in dir(fcts) if not d.startswith('_')]) 
7675              fun = None 
7676          elif argv[1] in ('install', '--install'): 
7677              from distutils.core import setup 
7678              setup(name="cma", 
7679                    long_description=__doc__, 
7680                    version=__version__.split()[0], 
7681                    description="CMA-ES, Covariance Matrix Adaptation Evolution Strategy for non-linear numerical optimization in Python", 
7682                    author="Nikolaus Hansen", 
7683                    author_email="hansen at lri.fr", 
7684                    maintainer="Nikolaus Hansen", 
7685                    maintainer_email="hansen at lri.fr", 
7686                    url="https://www.lri.fr/~hansen/cmaes_inmatlab.html#python", 
7687                    license="MIT", 
7688                    classifiers = [ 
7689                      "Intended Audience :: Science/Research", 
7690                      "Intended Audience :: Education", 
7691                      "Intended Audience :: Other Audience", 
7692                      "Topic :: Scientific/Engineering", 
7693                      "Topic :: Scientific/Engineering :: Mathematics", 
7694                      "Topic :: Scientific/Engineering :: Artificial Intelligence", 
7695                      "Operating System :: OS Independent", 
7696                      "Programming Language :: Python :: 2.6", 
7697                      "Programming Language :: Python :: 2.7", 
7698                      "Programming Language :: Python :: 3", 
7699                      "Development Status :: 4 - Beta", 
7700                      "Environment :: Console", 
7701                      "License :: OSI Approved :: MIT License", 
7702                    ], 
7703                    keywords=["optimization", "CMA-ES", "cmaes"], 
7704                    py_modules=["cma"], 
7705                    requires=["numpy"], 
7706              ) 
7707              fun = None 
7708          elif argv[1] in ('plot',): 
7709              plot(name=argv[2] if len(argv) > 2 else None) 
7710              raw_input('press return') 
7711              fun = None 
7712          elif len(argv) > 3: 
7713              fun = eval('fcts.' + argv[1]) 
7714          else: 
7715              print('try -h option') 
7716              fun = None 
7717   
7718      if fun is not None: 
7719   
7720          if len(argv) > 2:  # dimension 
7721              x0 = np.ones(eval(argv[2])) 
7722          if len(argv) > 3:  # sigma 
7723              sig0 = eval(argv[3]) 
7724   
7725          opts = {} 
7726          for i in xrange(5, len(argv), 2): 
7727              opts[argv[i - 1]] = eval(argv[i]) 
7728   
7729          # run fmin 
7730          if fun is not None: 
7731              tic = time.time() 
7732              fmin(fun, x0, sig0, opts)  # ftarget=1e-9, tolfacupx=1e9, verb_log=10) 
7733              # plot() 
7734              # print ' best function value ', res[2]['es'].best[1] 
7735              print('elapsed time [s]: + %.2f', round(time.time() - tic, 2)) 
7736   
7737      elif not len(argv): 
7738          fmin(fcts.elli, np.ones(6) * 0.1, 0.1, {'ftarget':1e-9}) 
7739   
7740   
7741  # ____________________________________________________________ 
7742  # ____________________________________________________________ 
7743  # 
7744  # mainly for testing purpose 
7745  # executed when called from an OS shell 
7746  if __name__ == "__main__": 
7747      # for i in range(1000):  # how to find the memory leak 
7748      #     main(["cma.py", "rastrigin", "10", "5", "popsize", "200", "maxfevals", "24999", "verb_log", "0"]) 
7749      main() 
7750   