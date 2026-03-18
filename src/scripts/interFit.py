#!/usr/bin/env python3

"""interFit.py -- Gaussian fitting of Grasp2K relativistic orbitals.

Reads relativistic orbital data from Grasp2K (via org_radwavefn), performs a
parameter sweep over maximum Gaussian exponent and number of terms per angular
momentum channel, and finds the expansion with the lowest weighted RMSE.

The fitting uses a geometric series of Gaussian exponents (alphas) and solves
the minimum-norm least squares system  A^T A x = A^T B  via pivotless LU
decomposition, where A[i,j] = exp(-alpha_j * r_i^2) and B is the matrix of
orbital function values on the radial grid.

Usage::

    interFit.py -comp 1              # single-component (averaged j)
    interFit.py -comp 2              # two-component (separate j)
    interFit.py -comp 4              # four-component (large + small)
    interFit.py -comp 1 -shrink 5 1  # with shrinking function (rc=5, sigma=1)
    interFit.py -tol 1e-4            # tighter fitting tolerance

Performance notes
-----------------
The following optimizations have been identified but not yet implemented:

1. **Vectorize get_A**: Replace the nested Python loop with numpy broadcasting:
   ``np.exp(-alphas[np.newaxis,:] * grid[:,np.newaxis]**2)``.
   The same applies to the RMSE summation in get_prefactor_coefficients.

2. **Vectorize pivotless_lu_decomp**:
   - Replace ``copy.deepcopy(A)`` with ``A.copy()`` (deepcopy pickles/unpickles
     numpy arrays, which is far slower than a direct memory copy).
   - Move the ``U[j,j]==0`` guard before the division (currently the division
     executes first, silently producing inf/nan).
   - Fold L into U in-place rather than allocating a separate L and adding at
     the end (L's diagonal is zero and the two triangles don't overlap).
   - Wrap the function with ``@numba.njit``.  The matrix is small (max 25x25)
     but the function is called ~360,000 times; Numba eliminates the Python
     loop overhead and typically gives 50-100x speedup for this pattern.
   - Note: scipy.linalg.lu_factor cannot be used as a replacement because it
     always applies partial pivoting (LAPACK dgetrf), which would disrupt the
     specific structure of the A^T A matrix in this application.

3. **Avoid unnecessary deep copies**: In get_prefactor_coefficients, the
   ``copy.deepcopy`` calls on numpy arrays (alphas, A, ATA, AT, functions)
   should use ``.copy()``.  The B matrix deep copy on the caller side
   (fitting_procedure) is also redundant since get_prefactor_coefficients
   already copies it.

4. **Parallelize the outer alpha loop**: The ``for maximum_alpha in
   self.loop_max_alpha`` loop is the natural parallelization target.  The
   itct early-exit counter resets per alpha value, so each alpha iteration is
   nearly independent.  The only cross-iteration dependency is the global
   best_RMSE/best_weight used for pruning; without it, workers do slightly
   more work but the result is correct.  Approach: split loop_max_alpha into
   N chunks (one per core), run each chunk with multiprocessing.Pool or
   concurrent.futures.ProcessPoolExecutor, then merge results by picking the
   global best across all workers.

Decoupling s and p term counts
------------------------------
Historically, s and p channels have always shared the same number of Gaussian
terms (the ``sp`` variable in fitting_procedure).  There is no physical reason
for this constraint, and allowing them to differ may improve fitting quality.
The constraint should be s >= p >= d >= f >= g.

Analysis of what would need to change:

- **Solver (get_prefactor_coefficients)**: Already l-aware.  ``num_terms[l]``
  is used per orbital (line ~564, ~576) and the LU solve block operates per-l
  with ``num_terms[i]``.  No changes needed.

- **Output (write_olcao_atom_type)**: Already writes s and p counts separately
  via ``num_terms[0]`` and ``num_terms[1]``.  No changes needed.

- **Sweep loop (fitting_procedure)**: Currently ``sp`` controls both s and p
  and drives the A matrix / LU construction.  With decoupled counts, p becomes
  an inner loop (p ranges from s down to min_num_terms).  The A matrix and LU
  decomposition depend only on s (the largest term count), so they do NOT need
  rebuilding when only p changes — only when s changes.

- **Weight function**: ``4*sp`` would split into separate s and p terms, e.g.
  ``2*s + 2*p`` (or another weighting to be determined).

- **_generate_term_sweeps**: p would become the outermost sweep channel
  (between s and d), maintaining s >= p >= d >= f >= g.

- **Search space**: Adding a separate p loop introduces up to 18 additional
  iterations per s value (~18x more work).  This makes the parallelization
  improvements (item 4 above) and the vectorization / Numba improvements
  (items 1-3) significantly more important.

- **RMSE loop concern**: In get_prefactor_coefficients (line ~587), the RMSE
  summation uses ``range(num_terms[0])`` (the s count) for ALL orbitals.  With
  s >= p this doesn't break, but the intent is unclear — it may be summing
  residuals over only the first s grid points rather than the full grid.  This
  should be verified against expected RMSE values before the s/p decoupling is
  implemented, as it may be a latent bug.
"""

import argparse as ap
import numpy as np
import scipy
import time
import copy
import os
import sys
import math
import pandas as pd
from datetime import datetime
try:
    import make_veusz_graph
except ImportError:
    make_veusz_graph = None
import org_radwavefn

from element_data import ElementData
from org_radwavefn import lktol, sktol, ltosk, ltok, ltolsym

from scipy.linalg import lu_solve

start=time.time()


# Define the main class that holds script data structures and settings.
class ScriptSettings():
    """The instance variables of this object are the user settings that
       control the program. The variable values are pulled from a list
       that is created within a resource control file and that are then
       reconciled with command line parameters."""


    def __init__(self):
        """Define default values for the parameters by pulling them
        from the resource control file in the default location:
        $OLCAO_RC/interFitrc.py or from the current working directory
        if a local copy of interFitrc.py is present."""

        # Read default variables from the resource control file.
        rc_dir = os.getenv('OLCAO_RC')
        if not rc_dir:
            sys.exit("Error: $OLCAO_RC is not set. See instructions.")
        sys.path.insert(1, rc_dir)
        from interFitrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults file.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the rc file.
        self.reconcile(args)

        # At this point, the command line parameters are set and accepted.
        #   When this initialization subroutine returns the script will
        #   start running. So, we use this as a good spot to record the
        #   command line parameters that were used.
        self.recordCLP()


    def assign_rc_defaults(self, default_rc):

        # Component mode default.
        self.number_components = default_rc["comp"]

        # Shrinking function defaults.
        self.shrink_rc = default_rc["shrink_rc"]
        self.shrink_sig = default_rc["shrink_sig"]
        self.tolerance = default_rc["tolerance"]


    def parse_command_line(self):

        # Create the parser tool.
        prog_name = "interFit.py"

        description_text = """
Gaussian fitting of Grasp2K relativistic orbitals for OLCAO.

Reads relativistic orbital data from Grasp2K (via org_radwavefn), performs a
parameter sweep over maximum Gaussian exponent and number of terms per angular
momentum channel, and finds the expansion with the lowest weighted RMSE.

The fitting uses a geometric series of Gaussian exponents (alphas) and solves
the minimum-norm least squares system  A^T A x = A^T B  via pivotless LU
decomposition, where A[i,j] = exp(-alpha_j * r_i^2) and B is the matrix of
orbital function values on the radial grid.

Operation modes:

  -comp 1/2/4
    Designates the amount of components to generate
     one, two, or four.
    1: Averages together to large components into a
       single orbital.
       -comp 1
    2: Fits only the large components.
       -comp 2
    3: Fits all four components.
       -comp 4

  -shrink rc sig
    Performs a shrinking function operation where rc
     is the critical radii and sig is the shrinking
     parameter.

           |          2
           |    (r-rc)
           |    ------
           |         2
           |    2*sig
     s(r)= |1-e         r .leq. rc
           |0           r > rc

     A file (shrink.inf) must be present with a list
      of all orbitals and corresponding rc sig in comma
      delimited columns.
"""

        epilog_text = """
Examples:
    interFit.py -comp 1              # single-component (averaged j)
    interFit.py -comp 2              # two-component (separate j)
    interFit.py -comp 4              # four-component (large + small)
    interFit.py -comp 1 -shrink 5 1  # with shrinking function (rc=5, sigma=1)
    interFit.py -tol 1e-4            # tighter fitting tolerance

Originally written by Ryan Thomas.
Please contact Paul Rulis (rulisp@umkc.edu) regarding questions.
Defaults are given in ./interFitrc.py or $OLCAO_RC/interFitrc.py.
"""

        parser = ap.ArgumentParser(prog = prog_name,
                formatter_class=ap.RawDescriptionHelpFormatter,
                description = description_text,
                epilog = epilog_text)

        # Add arguments to the parser.
        self.add_parser_arguments(parser)

        # Parse the arguments and return the results.
        return parser.parse_args()


    def add_parser_arguments(self, parser):

        # Define the component argument.
        parser.add_argument('-comp', dest='number_components',
                            type=int, default=self.number_components,
                            choices=[1, 2, 4],
                            help='Number of components to generate: '
                            '1=averaged j, 2=large components only, '
                            f'4=all four components. '
                            f'Default: {self.number_components}')

        # Define the shrinking function argument.
        parser.add_argument('-shrink', nargs=2, dest='shrink',
                            type=float, default=None,
                            metavar=('RC', 'SIG'),
                            help='Apply shrinking function with critical '
                            'radius RC and sigma SIG. '
                            f'Default: {self.shrink_rc} {self.shrink_sig}')

        # Define the tolerance argument.
        parser.add_argument('-tol', dest='tolerance',
                            type=float, default=self.tolerance,
                            metavar='TOL',
                            help='RMSE tolerance for fitting acceptance. '
                            f'Default: {self.tolerance}')


    def reconcile(self, args):
        self.number_components = args.number_components
        self.tolerance = args.tolerance

        if args.shrink is not None:
            self.shrink_rc = args.shrink[0]
            self.shrink_sig = args.shrink[1]
            self.apply_shrink = True
        else:
            self.apply_shrink = (self.shrink_rc != 0.0
                                 or self.shrink_sig != 0.0)


    def recordCLP(self):
        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime("%b. %d, %Y: %H:%M:%S")
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write(f"Cmnd:")
            for argument in sys.argv:
                cmd.write(f" {argument}")
            cmd.write("\n\n")


#==========================
# PRINTING FUNCTIONS:
#==========================

def print_line(number_lines):
  for i in range(number_lines):
    print("==================================================================================")

def print_callout(callout_text):
  print_line(1)
  print(callout_text)

def printStats(MaxA,numberTerms,Weight):
  print(str("%.8E"%MaxA) + "\t" + str(numberTerms[0]) + "\t" +
        str(numberTerms[1]) + "\t" + str(numberTerms[2]) + "\t" +
        str(numberTerms[3]) + "\t" + str(numberTerms[4]) + "\t" +
        str(Weight) + "\n")

def printProgress(minAlpha,maxAlpha,currentA,step):
  perc=((float(currentA)-float(minAlpha))/(float(maxAlpha)-float(minAlpha)))*100.0
  if perc%10==0:
    print(str("%.2f"%perc)+"%\t("+str(minAlpha)+"\t- "+str(currentA)+" -\t"+str(maxAlpha)+")")


#==========================
# ELEMENT INFO FUNCTIONS:
#==========================
def get_max_alpha(Z):
  """Return the maximum Gaussian exponent (alpha) for relativistic fitting.

  This table is specific to Grasp2K → OLCAO relativistic fitting and differs
  from the non-relativistic max_term_wf values in elements.dat.  Values
  increase with Z because heavier atoms have more tightly bound core orbitals
  requiring steeper Gaussians.
  """
  alphas={'1': 10000.0, '2': 10000.0, '3': 50000.0, '4': 50000.0, '5': 50000.0, '6': 50000.0, '7': 50000.0, '8': 50000.0, '9': 50000.0, '10': 50000.0, '11': 100000.0, '12': 100000.0, '13': 100000.0, '14': 200000.0, '15': 100000.0, '16': 100000.0, '17': 100000.0, '18': 100000.0, '19': 500000.0, '20': 500000.0, '21': 500000.0, '22': 500000.0, '23': 500000.0, '24': 500000.0, '25': 500000.0, '26': 500000.0, '27': 500000.0, '28': 500000.0, '29': 500000.0, '30': 1000000.0, '31': 1000000.0, '32': 1000000.0, '33': 1000000.0, '34': 1000000.0, '35': 1000000.0, '36': 1000000.0, '37': 5000000.0, '38': 5000000.0, '39': 5000000.0, '40': 5000000.0, '41': 5000000.0, '42': 5000000.0, '43': 5000000.0, '44': 5000000.0, '45': 5000000.0, '46': 5000000.0, '47': 5000000.0, '48': 5000000.0, '49': 5000000.0, '50': 5000000.0, '51': 5000000.0, '52': 5000000.0, '53': 5000000.0, '54': 5000000.0, '55': 10000000.0, '56': 10000000.0, '57': 10000000.0, '58': 10000000.0, '59': 10000000.0, '60': 10000000.0, '61': 10000000.0, '62': 10000000.0, '63': 10000000.0, '64': 10000000.0, '65': 10000000.0, '66': 10000000.0, '67': 10000000.0, '68': 10000000.0, '69': 10000000.0, '70': 10000000.0, '71': 50000000.0, '72': 50000000.0, '73': 50000000.0, '74': 50000000.0, '75': 50000000.0, '76': 50000000.0, '77': 50000000.0, '78': 50000000.0, '79': 50000000.0, '80': 50000000.0, '81': 50000000.0, '82': 50000000.0, '83': 50000000.0, '84': 50000000.0, '85': 50000000.0, '86': 50000000.0, '87': 50000000.0, '88': 50000000.0, '89': 50000000.0, '90': 50000000.0, '91': 50000000.0, '92': 50000000.0, '93': 50000000.0, '94': 50000000.0, '95': 50000000.0, '96': 50000000.0, '97': 50000000.0, '98': 50000000.0, '99': 50000000.0, '100': 50000000.0, '101': 50000000.0, '102': 50000000.0, '103': 50000000.0,}

  return alphas[str(Z)]


#==========================
# FORMATTING FUNCTIONS:
#==========================

def fortran_e(v):
  """Format a float in Fortran E18.8 style (e.g. '    0.12000000E+00').

  Produces the same output format as Fortran's E18.8 descriptor:
  a leading-zero mantissa (0.dddddddd) with an explicit exponent.
  The result is right-justified to 18 characters, matching the
  column layout used by the contract program's waveFn.dat writer.

  Parameters
  ----------
  v : float
      Value to format.

  Returns
  -------
  str
      18-character right-justified Fortran-style E-notation string.
  """
  if v == 0.0:
    return '0.00000000E+00'.rjust(18)
  sign = '-' if v < 0 else ''
  av = abs(v)
  exp = int(math.floor(math.log10(av))) + 1
  mantissa = av / (10.0 ** exp)
  mstr = f'{mantissa:.8f}'
  # Guard against rounding producing mantissa >= 1.0
  if float(mstr) >= 1.0:
    exp += 1
    mantissa /= 10.0
    mstr = f'{mantissa:.8f}'
  return f'{sign}{mstr}E{exp:+03d}'.rjust(18)


#==========================
# FITTING FUNCTIONS:
#==========================

def get_alpha_list(amin,amax,n,mode):
  """Generate a geometric series of Gaussian exponents from amin to amax.

  Parameters
  ----------
  amin : float
      Minimum exponent value.
  amax : float
      Maximum exponent value.
  n : int
      Number of terms in the series.
  mode : str
      Series type; currently only 'geometric' is supported.

  Returns
  -------
  list of float
      Exponent values: amin * (amax/amin)^((i-1)/(n-1)) for i=1..n.
  """
  a=[]
  if mode=='geometric':
    for i in range(1,n+1):
      r=amin*((float(amax/amin))**((i-1.0)/(n-1.0)))
      a.append(r)
    return a


def get_A(numS,radialGrid,GCoeffs):
  """Construct the A matrix for least squares: A[i,j] = exp(-alpha_j * r_i^2).

  Parameters
  ----------
  numS : int
      Number of Gaussian terms (columns).
  radialGrid : list of float
      Radial grid points (rows).
  GCoeffs : list of float
      Gaussian exponent values (alphas).

  Returns
  -------
  numpy.ndarray
      Matrix of shape (len(radialGrid), numS).
  """
  AMAT=np.zeros((len(radialGrid),numS),dtype="d")

  for j in range(numS):
    for i in range(len(radialGrid)):
      AMAT[i,j]=np.exp(-1.0*GCoeffs[j]*(radialGrid[i]**2.0))

  return AMAT


def pivotless_lu_decomp(A):
  """Pivotless LU decomposition returning L+U combined matrix.

  Performs Gaussian elimination without pivoting.  The result is a single
  matrix containing both L (lower triangle, unit diagonal implicit) and U
  (upper triangle), suitable for use with scipy.linalg.lu_solve when
  paired with a trivial pivot vector.

  Parameters
  ----------
  A : numpy.ndarray
      Square matrix to decompose.

  Returns
  -------
  numpy.ndarray
      Combined L+U matrix.
  """
  dim=np.shape(A)[0]

  L=np.zeros((dim,dim),dtype="d")
  U=copy.deepcopy(A)

  for j in range(dim-1):
    for i in range(j,dim-1):
      tc=(U[i+1,j]/U[j,j])
      if U[j,j]==0 or tc==0:
        continue
      else:
        U[i+1,:]=U[i+1,:]-tc*U[j,:]
        L[i+1,j]=tc

        # Adding this in to make truly upper triangular
        if U[i+1,j]<=1e-13 and U[i+1,j]>=-1e-13:
          U[i+1,j]=0.0

  return L+U


class orbital_fitting:
  """Gaussian fitting of Grasp2K relativistic orbitals to OLCAO basis functions.

  Performs a parameter sweep over the maximum Gaussian exponent and number of
  terms per angular momentum channel (s/p, d, f, g) to find the expansion
  that minimizes the total weighted RMSE.

  After fitting, results are stored as:

  - ``best_coefficients``: numpy matrix of fitted Gaussian coefficients
  - ``best_alphas``: list of Gaussian exponent values
  - ``best_num_terms``: [sp, sp, d, f, g] term counts
  - ``best_orbital_RMSEs``: per-orbital RMSE values
  - ``best_total_RMSE``: total weighted RMSE
  """

  def __init__(self, settings):
    self.apply_shrink=settings.apply_shrink
    self.shrink_critical_radius=settings.shrink_rc
    self.shrink_sigma=settings.shrink_sig
    self.number_components=settings.number_components

    self.element_data = ElementData()
    self.element_data.init_element_data()

    self.grasp_description=org_radwavefn.atomic_system('isodata', 'rwfn.out')

    self.min_alpha=0.12
    self.max_alpha=get_max_alpha(str(self.grasp_description.atomic_info.atomic_number))

    self.fitting_columns=['min_alpha', 'max_alpha', 'num_terms', 'num_s',
                          'num_p', 'num_d', 'num_f', 'num_g', 'weight',
                          'total_RMSE', 'RMSE', 'time']
    self.fitting_results=[]

    self.RMSE=1E100
    self.weight_factors=[1,3,5,7,9]
    self.weight=1E100
    self.fits=[]
    self.number_terms=[0,0,0,0,0]
    self.fitted_max_alpha=0

    self.orb_fit_tolerance=settings.tolerance
    self.alpha_step_size=int(self.max_alpha*0.001)

    self.loop_max_alpha=np.arange(0.3*self.max_alpha,
                                  20*self.max_alpha + self.alpha_step_size,
                                  self.alpha_step_size)
    self.max_num_terms=25
    self.min_num_terms=7

    # Best-fit results (populated by fitting_procedure)
    self.best_coefficients = None
    self.best_alphas = []
    self.best_num_terms = [0,0,0,0,0]
    self.best_orbital_RMSEs = []
    self.best_total_RMSE = 0.0


  def organize_orbitals(self):
    organize_time=time.time()

    if self.number_components==1:
      self.grasp_description.orbital_info.single_component()
    elif self.number_components==2:
      self.grasp_description.orbital_info.two_component()
    elif self.number_components==4:
      self.grasp_description.orbital_info.four_component()

    if self.apply_shrink:
      shrinking_time=time.time()
      self.grasp_description.orbital_info.apply_shrinking_function(self.shrink_critical_radius,self.shrink_sigma)
      print("  Shrinking Orbitals Complete:\t"+str('%.8E'%(time.time()-shrinking_time))+"s. / "+str('%.8E'%(time.time()-start))+"s.")
    else:
      for i in range(len(self.grasp_description.orbital_info.fitting_l_list)):
        for j in range(len(self.grasp_description.orbital_info.interpolated_radial_grid)):
          self.grasp_description.orbital_info.fitting_functions[i][j]=self.grasp_description.orbital_info.fitting_functions[i][j]/(self.grasp_description.orbital_info.interpolated_radial_grid[j]**float(self.grasp_description.orbital_info.fitting_l_list[i]))

    print("  Organizing Orbitals Complete:\t"+str('%.8E'%(time.time()-organize_time))+"s. / "+str('%.8E'%(time.time()-start))+"s.")


  def setup_calculation(self):
    calc_time=time.time()
    print(' Calculation Setup:')
    self.organize_orbitals()
    print(" Calculation Setup Complete:\t"+str('%.8E'%(time.time()-calc_time))+"s. / "+str('%.8E'%(time.time()-start))+"s.")


  def get_prefactor_coefficients(self,LU_matrix,A_matrix,AT_matrix,B_matrix,num_terms,alpha_list):
    """Solve the minimum-norm least squares system A^T A x = A^T B via LU.

    Each column of B is one orbital's function values on the radial grid.
    The system is solved block-wise: orbitals sharing the same l use the
    same number of Gaussian terms (and thus the same LU sub-block).

    Parameters
    ----------
    LU_matrix : numpy.ndarray
        LU decomposition of A^T A.
    A_matrix : numpy.ndarray
        The A matrix (Gaussians evaluated on the grid).
    AT_matrix : numpy.ndarray
        Transpose of A.
    B_matrix : numpy.ndarray
        Matrix of orbital function values (columns = orbitals).
    num_terms : list of int
        [sp, sp, d, f, g] number of terms per l-channel.
    alpha_list : list of float
        Gaussian exponent values.

    Returns
    -------
    tuple
        (coefficient_matrix, total_RMSE, per_orbital_RMSE, elapsed_time)
    """
    prefactor_time=time.time()

    alphas=copy.deepcopy(alpha_list)
    A=copy.deepcopy(A_matrix)
    ATA=copy.deepcopy(LU_matrix)
    AT=copy.deepcopy(AT_matrix)
    functions=copy.deepcopy(B_matrix)
    num_orbitals=copy.deepcopy(self.grasp_description.orbital_info.fitting_lvals)

    RMSE=[]

    function_number_terms=[]

    for l in self.grasp_description.orbital_info.fitting_l_list:
      function_number_terms.append(num_terms[l])

    AT_B=np.zeros((len(alphas),len(function_number_terms)),dtype="d")

    for j in range(len(function_number_terms)):
      for i in range(function_number_terms[j]):
        AT_B[i,j]=np.dot(AT[i,:],functions[:,j])

    coefficient_matrix=np.zeros((len(alphas),len(function_number_terms)),dtype="d")

    for i in range(len(num_orbitals)-1,-1,-1):
      if num_orbitals[i]!=0:
        x=num_terms[i]
        y1=sum(num_orbitals[i:])-num_orbitals[i]
        y2=y1+num_orbitals[i]
        PIV=np.arange(x)

        coefficient_matrix[:x,y1:y2]=scipy.linalg.lu_solve((ATA[:x,:x],PIV),AT_B[:x,y1:y2])

    RMSE=[0]*len(function_number_terms)
    RM=np.dot(A,coefficient_matrix)-functions

    for j in range(len(function_number_terms)):
      for i in range(num_terms[0]):
        RMSE[j]+=RM[i,j]**2.0
      RMSE[j]=np.sqrt(RMSE[j])

    total_RMSE=0.0

    for i in range(len(RMSE)):
      total_RMSE+=self.weight_factors[self.grasp_description.orbital_info.fitting_l_list[i]]*RMSE[i]

    return (coefficient_matrix,total_RMSE,RMSE,time.time()-prefactor_time)


  def _generate_term_sweeps(self, sp, fitting_max_l):
    """Generate all valid (d, f, g) term count combinations.

    Yields tuples of (d, f, g) where each is between sp (or previous
    channel) and min_num_terms, maintaining d >= f >= g.  For lower
    fitting_max_l, unused channels yield 0.

    Parameters
    ----------
    sp : int
        Number of terms for s and p channels.
    fitting_max_l : int
        Maximum angular momentum to fit (0=s, 1=p, 2=d, 3=f, 4=g).

    Yields
    ------
    tuple of (int, int, int)
        (d_terms, f_terms, g_terms)
    """
    if fitting_max_l <= 1:
      yield (0, 0, 0)
    elif fitting_max_l == 2:
      for d in range(sp, self.min_num_terms, -1):
        yield (d, 0, 0)
    elif fitting_max_l == 3:
      for d in range(sp, self.min_num_terms, -1):
        for f in range(d, self.min_num_terms, -1):
          yield (d, f, 0)
    elif fitting_max_l >= 4:
      for d in range(sp, self.min_num_terms, -1):
        for f in range(d, self.min_num_terms, -1):
          for g in range(f, self.min_num_terms, -1):
            yield (d, f, g)


  def fitting_procedure(self):
    """Sweep over max alpha and term counts to find optimal Gaussian expansion.

    Performs a parameter sweep over:
    1. Maximum alpha in the geometric series (0.3x to 20x of reference max_alpha)
    2. Number of terms for s/p channel (shared), then d, f, g channels

    The weight function penalizes more terms in higher-l channels:
    weight = 4*sp + 5*d + 7*f + 9*g.  Early-exit after 5 non-improving
    iterations per alpha value.

    Results are stored in self.best_coefficients, self.best_alphas,
    self.best_num_terms, self.best_orbital_RMSEs, and self.best_total_RMSE.
    """
    sweep_time=time.time()

    best_weight=self.weight
    best_RMSE=self.RMSE
    orbital_RMSEs=[]

    fitting_max_l = self.grasp_description.orbital_info.fitting_max_l

    # Channel display names
    channel_names = ['s', 'p', 'd', 'f', 'g']
    active_channels = channel_names[:fitting_max_l + 1]

    print_line(1)
    print(' Fitting Procedure:')
    print_line(1)
    print(f'Fitted Orbitals: {",".join(active_channels)}')
    print('Number of Orbitals:'
          + ''.join(f' {channel_names[l]}:{self.grasp_description.orbital_info.fitting_lvals[l]}'
                    for l in range(min(fitting_max_l + 1, 5))))

    best_coefficients = None
    best_alphas = []

    for maximum_alpha in self.loop_max_alpha:

      itct=0
      for sp in range(self.max_num_terms,self.min_num_terms,-1):

        fit_loop_time=time.time()

        if itct>=5:
          # For s/p-only (fitting_max_l<=1) the original code used 'break'.
          # For higher l the original used 'continue' on sp but 'break' on
          # inner loops.  Using 'break' here is equivalent since itct stays
          # >= 5 for all remaining sp values anyway.
          break

        alphas=get_alpha_list(self.min_alpha,maximum_alpha,sp,'geometric')
        A=get_A(sp, self.grasp_description.orbital_info.interpolated_radial_grid,alphas)
        AT=np.transpose(A)
        ATA=np.dot(AT,A)
        LU=pivotless_lu_decomp(ATA)
        B=np.transpose(np.array(copy.deepcopy(
            self.grasp_description.orbital_info.fitting_functions)))

        for d, f, g in self._generate_term_sweeps(sp, fitting_max_l):
          if itct>=5:
            break

          current_weight=4*sp+5*d+7*f+9*g

          if current_weight > best_weight:
            continue

          current_num_terms=[sp,sp,d,f,g]

          coefficients,total_RMSE,current_RMSE,fit_loop_time = \
                  self.get_prefactor_coefficients(LU,A,AT,B,current_num_terms,alphas)

          orbital_RMSEs.append(current_RMSE)

          if total_RMSE < best_RMSE or abs(total_RMSE - best_RMSE) <= self.orb_fit_tolerance:
            if round((maximum_alpha/max(self.loop_max_alpha))*100,2)%20==0:
              progress_args = [str('%.2f'%((maximum_alpha/max(self.loop_max_alpha))*100))+'%',
                    maximum_alpha, max(self.loop_max_alpha), sp]
              if fitting_max_l >= 2:
                progress_args.append(d)
              if fitting_max_l >= 3:
                progress_args.append(f)
              if fitting_max_l >= 4:
                progress_args.append(g)
              progress_args.append('%.5E'%total_RMSE)
              print(*progress_args)
            best_weight=current_weight
            best_RMSE=total_RMSE
            self.fitted_max_alpha=maximum_alpha
            self.number_terms=current_num_terms
            best_coefficients=coefficients
            best_alphas=copy.deepcopy(alphas)
            itct=0
          else:
            itct+=1

          self.fits.append([0.12,maximum_alpha,sp,sp,sp,d,f,g,
                            current_weight,total_RMSE,current_RMSE,
                            time.time()-fit_loop_time])

    # Store best results
    self.best_coefficients = best_coefficients
    self.best_alphas = best_alphas
    self.best_num_terms = self.number_terms
    self.best_orbital_RMSEs = orbital_RMSEs[-1] if orbital_RMSEs else []
    self.best_total_RMSE = best_RMSE

    # Write RMSE data
    if orbital_RMSEs:
      orb_frame_RMSE=pd.DataFrame(data=orbital_RMSEs,
          columns=self.grasp_description.orbital_info.fitting_orbital_names)
      orb_frame_RMSE.insert(0, 'iteration', range(1, len(orb_frame_RMSE) + 1))
      orb_frame_RMSE.to_csv('RMSE.plot', index=False, sep=' ')

    self.fitting_results=pd.DataFrame(data=self.fits,
                                      columns=self.fitting_columns)

    self.fitting_results.to_csv('fitting_results.plot', index=False, sep=' ')

    print(self.fitting_results.describe())

    if best_coefficients is not None:
      print("coefficients")
      print(best_coefficients)
      print("alphas")
      print(best_alphas)

    print(" Fitting Orbitals Complete:\t"+str('%.8E'%(time.time()-sweep_time))+"s. / "+str('%.8E'%(time.time()-start))+"s.")


  #==========================
  # OLCAO OUTPUT WRITER:
  #==========================

  def _parse_orbital_info(self):
    """Parse fitting results into structured orbital metadata.

    Returns a list of dicts sorted by (l ascending, n ascending) to match
    the olcao.dat orbital ordering used by the contract program.  Each
    dict has keys: 'name', 'n', 'l', 'coeff_index'.

    Returns
    -------
    list of dict
    """
    orbitals = []
    names = self.grasp_description.orbital_info.fitting_orbital_names
    l_list = self.grasp_description.orbital_info.fitting_l_list

    for idx, (name, l) in enumerate(zip(names, l_list)):
      # Names are like "1s", "2p", "3d", "10s" — everything up to the
      # last character is n, the last character is the l symbol.
      n = int(name[:-1])
      orbitals.append({'name': name, 'n': n, 'l': l, 'coeff_index': idx})

    orbitals.sort(key=lambda o: (o['l'], o['n']))
    return orbitals


  def _classify_orbitals(self, orbitals, Z):
    """Classify orbitals as core or valence and assign basis codes.

    Uses element_data core_orbitals[Z] to determine which orbitals are
    core (basisCode=1, present in all basis sets) and which are valence.
    Valence orbitals are assigned basis codes 1 (MB), 2 (FB), or 3 (EB)
    based on element_data vale_orbitals counts.

    Parameters
    ----------
    orbitals : list of dict
        Sorted orbital info from _parse_orbital_info.
    Z : int
        Atomic number.

    Returns
    -------
    tuple of (list of dict, list of dict)
        (core_orbitals, valence_orbitals), each with added 'basis_code'.
    """
    ed = self.element_data
    core_counts = ed.core_orbitals[Z]  # [n_s, n_p, n_d, n_f]

    l_seen = {}
    core_orbs = []
    vale_orbs = []

    for orb in orbitals:
      l = orb['l']
      l_seen[l] = l_seen.get(l, 0) + 1

      if l <= 3 and l_seen[l] <= core_counts[l]:
        orb['basis_code'] = 1
        core_orbs.append(orb)
      else:
        # Determine basis code from vale_orbitals counts
        vale_idx = l_seen[l] - (core_counts[l] if l <= 3 else 0)
        mb = ed.vale_orbitals[1][Z][l] if l <= ed.max_qn_l else 0
        fb = ed.vale_orbitals[2][Z][l] if l <= ed.max_qn_l else 0
        if vale_idx <= mb:
          orb['basis_code'] = 1
        elif vale_idx <= mb + fb:
          orb['basis_code'] = 2
        else:
          orb['basis_code'] = 3
        vale_orbs.append(orb)

    return core_orbs, vale_orbs


  def write_olcao_atom_type(self, f):
    """Write the OLCAO atom type section to a file object.

    Writes NUM_ATOM_TYPES, ATOM_TYPE_ID, ATOM_TYPE_LABEL, NUM_ALPHA_S_P_D_F,
    ALPHAS, core radial functions, and valence radial functions in the format
    expected by atomicTypes.f90.  Quantum number encoding follows
    contract/printResults.f90: QN_2j = 2*l (non-relativistic/single-component),
    numStates = (2l+1)*2, componentIndex = 1.

    Parameters
    ----------
    f : file-like object
        Open file to write to.
    """
    Z = self.grasp_description.atomic_info.atomic_number
    ed = self.element_data
    elem_name = ed.element_names[Z]
    num_terms = self.best_num_terms  # [sp, sp, d, f, g]
    alphas = self.best_alphas

    orbitals = self._parse_orbital_info()
    core_orbs, vale_orbs = self._classify_orbitals(orbitals, Z)

    # Header
    f.write("NUM_ATOM_TYPES\n")
    f.write("1    \n")
    f.write("ATOM_TYPE_ID__SEQUENTIAL_NUMBER\n")
    f.write("1     1     1         1\n")
    f.write("ATOM_TYPE_LABEL\n")
    f.write(f"{elem_name}1_1\n")

    # Alpha counts: s, p, d, f (drop g; lAngMomCount=4 in Fortran)
    f.write("NUM_ALPHA_S_P_D_F\n")
    f.write(f"{num_terms[0]:>8}{num_terms[1]:>8}{num_terms[2]:>8}{num_terms[3]:>8}\n")

    # Alphas (4 per line, Fortran E18.8 format)
    f.write("ALPHAS\n")
    for i in range(0, len(alphas), 4):
      chunk = alphas[i:i+4]
      f.write("".join(fortran_e(a) for a in chunk) + "\n")

    # Core radial functions
    nc = len(core_orbs)
    f.write("NUM_CORE_RADIAL_FNS\n")
    f.write(f"{nc:>8}{nc:>8}{nc:>8}\n")

    if nc > 0:
      f.write("NL_RADIAL_FUNCTIONS\n")
      for orb in core_orbs:
        n, l = orb['n'], orb['l']
        f.write(f"{1:>8}{orb['basis_code']:>8}\n")
        f.write(f"{n:>8}{l:>8}{2*l:>8}{(2*l+1)*2:>8}{1:>8}\n")
        nc_l = num_terms[l]
        idx = orb['coeff_index']
        coeffs = [float(self.best_coefficients[j, idx]) for j in range(nc_l)]
        for ci in range(0, len(coeffs), 4):
          chunk = coeffs[ci:ci+4]
          f.write("".join(fortran_e(c) for c in chunk) + "\n")

    # Valence radial functions (cumulative counts per basis level)
    mb_vale = sum(1 for o in vale_orbs if o['basis_code'] <= 1)
    fb_vale = sum(1 for o in vale_orbs if o['basis_code'] <= 2)
    eb_vale = len(vale_orbs)
    f.write("NUM_VALE_RADIAL_FNS\n")
    f.write(f"{mb_vale:>8}{fb_vale:>8}{eb_vale:>8}\n")

    if eb_vale > 0:
      f.write("NL_RADIAL_FUNCTIONS\n")
      for orb in vale_orbs:
        n, l = orb['n'], orb['l']
        f.write(f"{1:>8}{orb['basis_code']:>8}\n")
        f.write(f"{n:>8}{l:>8}{2*l:>8}{(2*l+1)*2:>8}{1:>8}\n")
        nc_l = num_terms[l]
        idx = orb['coeff_index']
        coeffs = [float(self.best_coefficients[j, idx]) for j in range(nc_l)]
        for ci in range(0, len(coeffs), 4):
          chunk = coeffs[ci:ci+4]
          f.write("".join(fortran_e(c) for c in chunk) + "\n")


  def write_olcao_pot_type(self, f):
    """Write the OLCAO potential type section to a file object.

    Writes NUM_POTENTIAL_TYPES, POTENTIAL_TYPE_ID, POTENTIAL_TYPE_LABEL,
    NUCLEAR_CHARGE__ALPHA, COVALENT_RADIUS, NUM_ALPHAS, and ALPHAS in the
    format expected by potTypes.f90.  The nuclear alpha is always 20.0
    (standard OLCAO convention).  Potential alpha range and count come from
    element_data (num_terms_pot, min_term_pot, max_term_pot).

    Parameters
    ----------
    f : file-like object
        Open file to write to.
    """
    Z = self.grasp_description.atomic_info.atomic_number
    ed = self.element_data
    elem_name = ed.element_names[Z]

    f.write("NUM_POTENTIAL_TYPES\n")
    f.write("1    \n")
    f.write("POTENTIAL_TYPE_ID__SEQUENTIAL_NUMBER\n")
    f.write("1     1     1         1\n")
    f.write("POTENTIAL_TYPE_LABEL\n")
    f.write(f"{elem_name}1_1\n")
    f.write("NUCLEAR_CHARGE__ALPHA\n")
    f.write(f"{float(Z):f} {20.0:f}\n")
    f.write("COVALENT_RADIUS\n")
    f.write(f"{ed.coval_radii[Z]:f}\n")
    f.write("NUM_ALPHAS\n")
    f.write(f"{ed.num_terms_pot[Z]}\n")
    f.write("ALPHAS\n")
    f.write(f"{ed.min_term_pot[Z]:e} {ed.max_term_pot[Z]:e}\n")


  def write_olcao_snippet(self, filename=None):
    """Write combined atom type + potential type sections.

    This is the primary output method: it produces a file containing the
    atomic basis function and potential type data in the format read by
    OLCAO's atomicTypes.f90 and potTypes.f90.

    Parameters
    ----------
    filename : str, optional
        Output file path.  If None, writes to stdout.

    Raises
    ------
    RuntimeError
        If called before fitting_procedure has completed.
    """
    if self.best_coefficients is None:
      raise RuntimeError(
          "No fitting results available. Run fitting_procedure first.")

    if filename:
      with open(filename, 'w') as f:
        self.write_olcao_atom_type(f)
        self.write_olcao_pot_type(f)
      print(f"  OLCAO snippet written to {filename}")
    else:
      self.write_olcao_atom_type(sys.stdout)
      self.write_olcao_pot_type(sys.stdout)


#==========================
# PROGRAM OPERATIONS:
#==========================

def main():

    # Get script settings from a combination of the resource control file
    #   and parameters given by the user on the command line.
    settings = ScriptSettings()

    print_line(1)
    print_callout("\t\tI N T E R F I T    F I T T I N G    P R O G R A M")
    print_line(1)

    # Start executing the main activities of the program.
    element_system=orbital_fitting(settings)

    element_system.setup_calculation()

    element_system.fitting_procedure()

    element_system.write_olcao_snippet('waveFn.dat')

    # Finalize the program activities and quit.
    print_callout("Total Time: "+str('%.8E'%(time.time()-start))+"s.")
    print_line(1)
    print_line(1)


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or a request
    #   to import information from external modules. Only now do we actually
    #   start running the program. The purpose of this is to allow another
    #   python program to import *this* script and call its functions
    #   internally.
    main()
