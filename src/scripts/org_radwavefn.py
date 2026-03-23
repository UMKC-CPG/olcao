"""Parser for Grasp2K Dirac 4-component relativistic wavefunctions (rwfn.out).

Reads Grasp2K output (rwfn.out) containing numerical radial wavefunctions
P(r) (large component) and Q(r) (small component) on a log-linear grid,
and organizes them for Gaussian fitting.  Also parses the companion
``isodata`` file for nuclear charge-distribution parameters.

The module provides three fitting modes that trade off relativistic
detail for simplicity:

- **single_component**: averages j-split pairs into one orbital per (n,l)
- **two_component**: keeps each j-component separate
- **four_component**: fits both large and small Dirac components per j

A Gaussian damping / shrinking function can be applied before fitting to
ensure the orbitals vanish at a finite cutoff radius.
"""

import time
import numpy as np
import math
import copy
import pandas as pd
import matplotlib.pyplot as plt


# Dictionaries for Z to element symbol and back.
el_num_dict={'H': 1 , 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb':37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd':46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs':55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109}
el_syb_dict={'1': 'H', '2': 'He', '3': 'Li', '4': 'Be', '5': 'B', '6': 'C', '7': 'N', '8': 'O', '9': 'F', '10': 'Ne', '11': 'Na', '12': 'Mg', '13': 'Al', '14': 'Si', '15': 'P', '16': 'S', '17': 'Cl', '18': 'Ar', '19': 'K', '20': 'Ca', '21': 'Sc', '22': 'Ti', '23': 'V', '24': 'Cr', '25': 'Mn', '26': 'Fe', '27': 'Co', '28': 'Ni', '29': 'Cu', '30': 'Zn', '31': 'Ga', '32': 'Ge', '33': 'As', '34': 'Se', '35': 'Br', '36': 'Kr', '37': 'Rb', '38': 'Sr', '39': 'Y', '40': 'Zr', '41': 'Nb', '42': 'Mo', '43': 'Tc', '44': 'Ru', '45': 'Rh', '46': 'Pd', '47': 'Ag', '48': 'Cd', '49': 'In', '50': 'Sn', '51': 'Sb', '52': 'Te', '53': 'I', '54': 'Xe', '55': 'Cs', '56': 'Ba', '57': 'La', '58': 'Ce', '59': 'Pr', '60': 'Nd', '61': 'Pm', '62': 'Sm', '63': 'Eu', '64': 'Gd', '65': 'Tb', '66': 'Dy', '67': 'Ho', '68': 'Er', '69': 'Tm', '70': 'Yb', '71': 'Lu', '72': 'Hf', '73': 'Ta', '74': 'W', '75': 'Re', '76': 'Os', '77': 'Ir', '78': 'Pt', '79': 'Au', '80': 'Hg', '81': 'Tl', '82': 'Pb', '83': 'Bi', '84': 'Po', '85': 'At', '86': 'Rn', '87': 'Fr', '88': 'Ra', '89': 'Ac', '90': 'Th', '91': 'Pa', '92': 'U', '93': 'Np', '94': 'Pu', '95': 'Am', '96': 'Cm', '97': 'Bk', '98': 'Cf', '99': 'Es', '100': 'Fm', '101': 'Md', '102': 'No', '103': 'Lr', '104': 'Rf', '105': 'Db', '106': 'Sg', '107': 'Bh', '108': 'Hs', '109': 'Mt'}


# Functions to convert between the Dirac kappa quantum number and
# angular momentum l. In the Dirac equation, each non-relativistic
# orbital with angular momentum l splits into two states with total
# angular momentum j = l +/- 1/2. The kappa quantum number encodes
# both l and j:
#   kappa = -(l+1)  for j = l + 1/2
#   kappa =  l      for j = l - 1/2
# The large and small components of the Dirac spinor have different
# angular momentum: if the large component has l, the small component
# has l' = l +/- 1 (depending on the sign of kappa).


def lktol(k):
  """Convert Dirac kappa to angular momentum l for the large component.

  For kappa < 0 (j = l + 1/2): l = |kappa| - 1
  For kappa > 0 (j = l - 1/2): l = kappa

  Examples: kappa=-1 -> l=0 (s), kappa=1 -> l=1 (p1/2),
            kappa=-2 -> l=1 (p3/2), kappa=2 -> l=2 (d3/2)
  """
  if k<0:
    return (-1*k)-1
  else:
    return k

def sktol(k):
  """Convert Dirac kappa to angular momentum l for the small component.

  The small component of the Dirac spinor has angular momentum that
  differs from the large component by one unit:
    For kappa < 0: l_small = |kappa|     (= l_large + 1)
    For kappa > 0: l_small = kappa - 1   (= l_large - 1)

  Examples: kappa=-1 -> l_small=1 (small comp of 1s is p-like),
            kappa=1 -> l_small=0 (small comp of 2p1/2 is s-like)
  """
  if k<0:
    return k*-1
  else:
    return k-1

def ltosk(l):
  """Convert angular momentum l to both kappa values for the small component.

  Returns the two kappa values whose small component has angular
  momentum l. For l=0, only kappa=-1 has an s-like small component,
  so (1, 1) is returned (both map to the same state). For l>0,
  returns (kappa1, kappa2) = (-l, l+1).

  Returns:
      tuple: (kappa1, kappa2) where sktol(kappa_i) == l.
  """
  if l==0:
    return 1,1
  else:
    return -1*l,l+1

def ltok(l):
  """Convert angular momentum l to both kappa values for the large component.

  Returns the two kappa values whose large component has angular
  momentum l. For l=0 (s orbitals), only kappa=-1 exists (j=1/2),
  so (-1, -1) is returned. For l>0, returns (kappa1, kappa2):
    kappa1 = l       (j = l - 1/2)
    kappa2 = -(l+1)  (j = l + 1/2)

  Returns:
      tuple: (kappa1, kappa2) where lktol(kappa_i) == l.
  """
  if l==0:
    return -1,-1
  else:
    return l,-1*(l+1)

def ltolsym(l):
  """Convert angular momentum quantum number l to spectroscopic symbol.

  Mapping: 0->'s', 1->'p', 2->'d', 3->'f', 4->'g'
  """
  if l==0:
    return "s"
  elif l==1:
    return "p"
  elif l==2:
    return "d"
  elif l==3:
    return "f"
  elif l==4:
    return "g"


class isodata_file:
  """Nuclear charge distribution parameters from a Grasp2K isodata file.

  Parses the ``isodata`` file produced by Grasp2K, which describes the
  Fermi-model nuclear charge distribution.  Key parameters:

  - **fermi_dist_param_a**: diffuseness parameter *a* (fm) of the Fermi
    distribution rho(r) = rho0 / (1 + exp((r - c) / a)).
  - **fermi_dist_param_c**: half-density radius *c* (fm), the radius at
    which the nuclear charge density drops to half its central value.

  Also records: atomic number *Z*, mass number *A*, nuclear mass (amu),
  nuclear spin *I* (units of hbar), magnetic dipole moment (nuclear
  magnetons), and electric quadrupole moment (barns).
  """

  def parse_file(self,fname):
    with open(fname,'r') as f:
      file=f.readlines()
    file=[x.strip().split() for x in file]
    
    self.atomic_number=int(float(file[1][0]))
    self.mass_number=int(float(file[3][0]))
    self.fermi_dist_param_a=float(file[5][0])
    self.fermi_dist_param_c=float(file[7][0])
    self.nucleus_mass=float(file[9][0])
    self.nuclear_spin=float(file[11][0])
    self.nuclear_dipole_moment=float(file[13][0])
    self.nuclear_quadrupole_moment=float(file[15][0])


  def __init__(self,fname):
    self.atomic_number=0
    self.mass_number=0  
    self.fermi_dist_param_a=0.0
    self.fermi_dist_param_c=0.0
    self.nucleus_mass=0.0 #in amu
    self.nuclear_spin=0.0   #(I) (in units of h / 2 pi)
    self.nuclear_dipole_moment=0.0 #(in nuclear magnetons)
    self.nuclear_quadrupole_moment=0.0 #(in barns)

    self.parse_file(fname)

class orbital_file:
  """Relativistic orbital wavefunctions parsed from a Grasp2K rwfn.out file.

  Reads the numerical radial wavefunctions (large component P(r) and small
  component Q(r)) for every orbital on the Grasp2K log-linear grid.  After
  parsing, provides interpolated data on a unified grid suitable for
  Gaussian fitting.
  """

  def parse_file(self,fname):
    """Parse a Grasp2K ``rwfn.out`` file and populate orbital data.

    File format:
      - Header line: ``"G92RWF"`` (ignored).
      - For each orbital: a 4-field header ``n  kappa  eigenvalue  num_points``,
        followed by a ZORA coefficient on the next line, then three blocks of
        *num_points* values: P(r) large component, Q(r) small component, and
        the radial grid points r.
      - Fortran ``'D'`` exponents are converted to ``'E'`` for Python parsing.
      - The first grid point (r = 0) is dropped because subsequent processing
        divides by r.

    After parsing, interpolated arrays on a merged radial grid and
    component-divided-by-r arrays are constructed for all orbitals.
    """

    # Open the file and store lines as a list.
    with open(fname,'r') as f:
      rwfn=f.readlines()

      # FORTRAN to PYTHON - Scientific notation should be 'E'
      for i in range(len(rwfn)):
      	rwfn[i]=rwfn[i].replace('D','E').strip()

    # Loop through the rwfn file and when each line has four items in it
    #   then an orbital description is starting. What follows is three
    #   segments of large components, small components, and radial
    #   components. After these are parsed and added to the respective
    #   increment the loop variable.
    i=0
    while i <len(rwfn):
      
      l=rwfn[i].split()

      # Begin an orbital description
      if len(l)==4:
        self.num_orbitals+=1
        self.n_qnums.append(int(l[0]))
        self.k_qnums.append(int(l[1]))
        self.qnums.append([int(l[0]),int(l[1]),lktol(int(l[1])),sktol(int(l[1]))])
        self.energy_eigenvalues.append(float(l[2]))
        self.numerical_points.append(int(l[3]))

        i+=1

        self.zero_point_energies.append(float(rwfn[i].replace('D','E').strip()))

        i+=1
        # Because there are three segments, we want to find the bounds of
        #  those lists. Also, each one starts with a 0.000000 which is
        #  problematic for dividing by r so we jettison that value here
        #  and start on the next non-zero value.
        n=self.numerical_points[-1]

        bounds=[i+1,(n+i)-1,(n+i)+1,2*n+i-1,2*n+i+1,3*n+i-1]

        self.large_comps.append([float(x) for x in rwfn[bounds[0]:bounds[1]]])
        self.small_comps.append([float(x) for x in rwfn[bounds[2]:bounds[3]]])
        self.radial_comps.append([float(x) for x in rwfn[bounds[4]:bounds[5]]])
        self.numerical_points[-1]=self.numerical_points[-1]-1
        

      else:
      	i+=1


    # Get maximum quantum numbers for fitting algorithms.
    self.max_n=max(self.n_qnums)
    self.max_k=max(self.k_qnums,key=abs)
    self.max_l=lktol(self.max_k)
    self.max_small_l=sktol(self.max_k)
    self.orbital_list=[str(x)+'_'+str(y) for x,y in zip(self.n_qnums,self.k_qnums)]

    # Make an interpolated radial grid for fitting purposes.
    self.interpolated_radial_grid=sorted(list(set([float(r) for sublist in self.radial_comps for r in sublist])))

    # Make the columns needed for the orbitals - positive and negative divided by r.   
    for i in range(self.num_orbitals):
      self.large_div_r.append([lc/rc for lc,rc in
          zip(self.large_comps[i],self.radial_comps[i])])
      self.neg_large_div_r.append([-1*lc/rc for lc,rc in
          zip(self.large_comps[i],self.radial_comps[i])])
      self.small_div_r.append([sc/rc for sc,rc in
          zip(self.small_comps[i],self.radial_comps[i])])
      self.neg_small_div_r.append([-1*sc/rc for sc,rc in
          zip(self.small_comps[i],self.radial_comps[i])])

      # Large component fitting.
      self.interpolated_large_components.append(
          np.interp(np.asarray(self.interpolated_radial_grid),
              np.asarray(self.radial_comps[i]),
              np.asarray(self.large_comps[i])).tolist())
      self.interpolated_ldivr_components.append([lc/rc for lc,rc in
          zip(self.interpolated_large_components[i],
              self.interpolated_radial_grid)])
      self.interpolated_nldivr_components.append([-1*lcr for lcr in
          self.interpolated_ldivr_components[i]])

      # Small component fitting.
      self.interpolated_small_components.append(
          np.interp(np.asarray(self.interpolated_radial_grid),
              np.asarray(self.radial_comps[i]),
              np.asarray(self.small_comps[i])).tolist())
      self.interpolated_sdivr_components.append([sc/rc for sc,rc in
          zip(self.interpolated_small_components[i],
              self.interpolated_radial_grid)])
      self.interpolated_nsdivr_components.append([-1*scr for scr in
          self.interpolated_sdivr_components[i]])

  def single_component(self):
    """Scalar-relativistic approximation: average j=l+1/2 and j=l-1/2
    large components into one non-relativistic-like orbital per (n,l).

    For l=0 (s orbitals), there is only one j=1/2 component so no
    averaging is needed.  Populates fitting_functions, fitting_lvals,
    fitting_l_list, and fitting_orbital_names.
    """
    self.fitting_max_l=self.max_l
    for j in range(self.max_l,-1,-1):
      for i in range(self.max_n,0,-1):
        if ltok(j)[0]==ltok(j)[1]:
          if str(i)+'_'+str(ltok(j)[0]) in self.orbital_list:
            self.fitting_orbital_names.append(str(i)+str(ltolsym(j)))
            self.fitting_functions.append(
                self.interpolated_large_components[
                    self.orbital_list.index(str(i)+'_'+str(ltok(j)[0]))])
            self.fitting_lvals[j]+=1
            self.fitting_l_list.append(j)
        else:
          if str(i)+'_'+str(ltok(j)[0]) in self.orbital_list \
              and str(i)+'_'+str(ltok(j)[1]) in self.orbital_list:
            self.fitting_orbital_names.append(str(i)+str(ltolsym(j)))
            self.fitting_functions.append([(f1+f2)/2.0 for f1,f2 in
                zip(self.interpolated_large_components[
                    self.orbital_list.index(str(i)+'_'+str(ltok(j)[0]))],
                    self.interpolated_large_components[
                        self.orbital_list.index(str(i)+'_'+str(ltok(j)[1]))])])
            self.fitting_lvals[j]+=1
            self.fitting_l_list.append(j)


  
  def two_component(self):
    """Two-component fitting: keeps each j-component separate.

    For example, 2p_1/2 and 2p_3/2 are fitted individually, preserving
    spin-orbit splitting information.  Each j-component's large-component
    P(r) is added to the fitting list as a separate function.
    """
    self.fitting_max_l=self.max_l
    for j in range(self.max_l,-1,-1):
      for i in range(self.max_n,0,-1):
        if ltok(j)[0]==ltok(j)[1]:
          if str(i)+'_'+str(ltok(j)[0]) in self.orbital_list:
            self.fitting_orbital_names.append(str(i)+str(ltolsym(j)))
            self.fitting_functions.append(self.interpolated_large_components[
                self.orbital_list.index(str(i)+'_'+str(ltok(j)[0]))])
            self.fitting_lvals[j]+=1
            self.fitting_l_list.append(j)
        else:
          if str(i)+'_'+str(ltok(j)[0]) in self.orbital_list:
            self.fitting_orbital_names.append(str(i)+str(ltolsym(j)))
            self.fitting_functions.append(self.interpolated_large_components[
                self.orbital_list.index(str(i)+'_'+str(ltok(j)[0]))])
            self.fitting_lvals[j]+=1
            self.fitting_l_list.append(j)

          if str(i)+'_'+str(ltok(j)[1]) in self.orbital_list:
            self.fitting_orbital_names.append(str(i)+str(ltolsym(j)))
            self.fitting_functions.append(self.interpolated_large_components[
                self.orbital_list.index(str(i)+'_'+str(ltok(j)[1]))])
            self.fitting_lvals[j]+=1
            self.fitting_l_list.append(j)

  def four_component(self):
    """Full Dirac 4-component: fits both P(r) large and Q(r) small
    components for each j, giving a complete relativistic description.

    Small components are prefixed with 's_' and large with 'l_' in the
    fitting orbital names.  The fitting_max_l is set to the maximum
    small-component angular momentum.
    """
    self.fitting_max_l=self.max_small_l
    for j in range(self.max_small_l,-1,-1):
      for i in range(self.max_n,0,-1):

        # Do small components first where now small component
        #  designation will have s_n_k.
        if ltosk(j)[0]==ltosk(j)[1]:
          if str(i)+'_'+str(ltosk(j)[0]) in self.orbital_list:
            self.fitting_orbital_names.append('s_'+str(i)+'_'+str(ltosk(j)[0]))
            self.fitting_functions.append(self.interpolated_small_components[
                self.orbital_list.index(str(i)+'_'+str(ltosk(j)[0]))])
            self.fitting_lvals[j]+=1
            self.fitting_l_list.append(j)
        else:
          if str(i)+'_'+str(ltosk(j)[0]) in self.orbital_list:
            self.fitting_orbital_names.append('s_'+str(i)+'_'+str(ltosk(j)[0]))
            self.fitting_functions.append(self.interpolated_small_components[
                self.orbital_list.index(str(i)+'_'+str(ltosk(j)[0]))])
            self.fitting_lvals[j]+=1
            self.fitting_l_list.append(j)
          if str(i)+'_'+str(ltosk(j)[1]) in self.orbital_list:
            self.fitting_orbital_names.append('s_'+str(i)+'_'+str(ltosk(j)[1]))
            self.fitting_functions.append(self.interpolated_small_components[
                self.orbital_list.index(str(i)+'_'+str(ltosk(j)[1]))])
            self.fitting_lvals[j]+=1
            self.fitting_l_list.append(j)

        # Do large components where now large components
        #   designation will be l_n_k.
        if ltok(j)[0]==ltok(j)[1]:
          if str(i)+'_'+str(ltok(j)[0]) in self.orbital_list:
            self.fitting_orbital_names.append('l_'+str(i)+'_'+str(ltok(j)[0]))
            self.fitting_functions.append(self.interpolated_large_components[
                self.orbital_list.index(str(i)+'_'+str(ltok(j)[0]))])
            self.fitting_lvals[j]+=1
            self.fitting_l_list.append(j)
        else:
          if str(i)+'_'+str(ltok(j)[0]) in self.orbital_list:
            self.fitting_orbital_names.append('l_'+str(i)+'_'+str(ltok(j)[0]))
            self.fitting_functions.append(self.interpolated_large_components[
                self.orbital_list.index(str(i)+'_'+str(ltok(j)[0]))])
            self.fitting_lvals[j]+=1
            self.fitting_l_list.append(j)
          if str(i)+'_'+str(ltok(j)[1]) in self.orbital_list:
            self.fitting_orbital_names.append('l_'+str(i)+'_'+str(ltok(j)[1]))
            self.fitting_functions.append(self.interpolated_large_components[
                self.orbital_list.index(str(i)+'_'+str(ltok(j)[1]))])
            self.fitting_lvals[j]+=1
            self.fitting_l_list.append(j)

  def apply_shrinking_function(self,rc,sigma):
    """Apply Gaussian damping to ensure orbitals vanish at finite radius.

    For r > rc, the orbital value is set to zero.  For r <= rc, the
    orbital is multiplied by a damping factor (1 - exp(-(r-rc)^2 / (2*sigma^2))).
    The orbital is also divided by r^l to extract the radial part for
    Gaussian fitting.

    Parameters
    ----------
    rc : float
        Critical radius (Bohr) beyond which the orbital is zeroed.
    sigma : float
        Width parameter for the Gaussian damping envelope.
    """
    for i in range(len(self.fitting_orbital_names)):
      for j in range(len(self.fitting_functions[i])):
        if self.interpolated_radial_grid[j] > rc:
          self.fitting_functions[i][j]=0.0
        else:
          self.fitting_functions[i][j] = \
              self.fitting_functions[i][j]*(
                  1-math.exp((-0.5*(self.interpolated_radial_grid[j]-rc)**2.0) /
                  (2.0*sigma**2.0)))
          self.fitting_functions[i][j] = \
              self.fitting_functions[i][j] / \
              self.interpolated_radial_grid[j]**np.float64(self.fitting_l_list[i])


  def __init__(self,fname):
    self.num_orbitals=0         # Number of orbitals in the system.
    self.n_qnums=[]             # List of principal quantum numbers.
    self.k_qnums=[]             # List of k quantum numbers.
    self.qnums=[]               # List of lists of quantum numbers [n,k,l_l,l_s].
    self.energy_eigenvalues=[]  # List of the energy eigenvalue for each orbital.
    self.numerical_points=[]    # List of the number of numerical points in each orbital definition.
    self.zero_point_energies=[] # List of the zero point energies for each orbital.
    self.large_comps=[]         # List of lists of large components.
    self.small_comps=[]         # List of lists of small components.
    self.radial_comps=[]        # List of lists of radial components.
    self.large_div_r=[]         # List of lists of large / radial components.
    self.neg_large_div_r=[]     # List of lists of -1 * large / radial components.
    self.small_div_r=[]         # List of lists of small / radial components.
    self.neg_small_div_r=[]     # List of lists of -1 * small / radial components.

    self.max_n=0                # Maximum of principal quantum number.
    self.max_l=0                # Maximum angular momentum quantum number.
    self.max_small_l=0          # Maximum angular momentum quantum number for small component.
    self.max_k=0                # Maximum of kappa quantum number.

    self.orbital_list=[]        # List to contain names of orbitals in n_k format.

    self.interpolated_radial_grid=[]
    self.interpolated_large_components=[]
    self.interpolated_small_components=[]
    self.interpolated_ldivr_components=[]
    self.interpolated_sdivr_components=[]
    self.interpolated_nldivr_components=[]
    self.interpolated_nsdivr_components=[]

    self.fitting_max_l=0
    self.fitting_orbital_names=[]
    self.fitting_functions=[]
    self.fitting_lvals=[0,0,0,0,0]
    self.fitting_l_list=[]

    self.parse_file(fname)


# Class that has the element symbol and automatically parses isodata
#   and rwfn.out files for contents.
class atomic_system:
  """Convenience class combining isodata and rwfn.out for a single element.

  Attributes
  ----------
  orbital_info : orbital_file
      Parsed relativistic orbital wavefunctions.
  atomic_info : isodata_file
      Nuclear charge distribution parameters.
  element : str
      Chemical symbol (e.g. 'C', 'Si').
  """

  def __init__(self,isofile,wavefile):
    self.orbital_info=orbital_file(wavefile)
    self.atomic_info=isodata_file(isofile)
    self.element=el_syb_dict[str(self.atomic_info.atomic_number)]


# If running the python script, create the classes and write out the desired
#  files as input to fitting and graphing solutions.
if __name__== "__main__":

  stime=time.time()

  # Print message that module file is running.
  print('Running Module File')

  # Create the atomic system class which will parse and calculate
  #  the required information about the atom from isodata and the
  #  orbital information from rwfn.out.
  atom=atomic_system('isodata','rwfn.out')

  # Create the eigenvalues file as a csv
  with open('eigenvalues','w') as f:
    f.write('n_quantum_number,k_quantum_number,large_l_quantum_number,small_l_quantum_number,energy_eigenvalue\n')    
    for i in range(atom.orbital_info.num_orbitals):
      f.write(str(atom.orbital_info.qnums[i][0]) + ','
          + str(atom.orbital_info.qnums[i][1]) + ','
          + str(atom.orbital_info.qnums[i][2]) + ','
          + str(atom.orbital_info.qnums[i][3]) + ','
          + str(atom.orbital_info.energy_eigenvalues[i]) + '\n')
  
  # Create a file for each orbital using <Element>_n_<nqnum>_k_<qnum>.dat 
  for i in range(atom.orbital_info.num_orbitals):
    with open(atom.element
              + '_n_' + str(atom.orbital_info.n_qnums[i])
              + '_k_' + str(atom.orbital_info.k_qnums[i])
              + '.dat', 'w') as g:

      g.write('radial_grid,large_component,large_component_div_r,'
          + 'neg_large_component_div_r,small_component,'
          + 'small_component_div_r,neg_small_component_div_r\n')
      for r, l, l_r, n_l_r, s, s_r, n_s_r in zip(
          atom.orbital_info.radial_comps[i],
          atom.orbital_info.large_comps[i],
          atom.orbital_info.large_div_r[i],
          atom.orbital_info.neg_large_div_r[i],
          atom.orbital_info.small_comps[i],
          atom.orbital_info.small_div_r[i],
          atom.orbital_info.neg_small_div_r[i]):
        g.write(str(r) + ','
            + str(l) + ','
            + str(l_r) + ','
            + str(n_l_r) + ','
            + str(s) + ','
            + str(s_r) + ','
            + str(n_s_r) + '\n')



  data_dict={}

  interpolated_data=[[copy.deepcopy(atom.orbital_info.interpolated_radial_grid)]]
  col_vals=['r']

  data_dict['r']=copy.deepcopy(atom.orbital_info.interpolated_radial_grid)

  for i in range(atom.orbital_info.num_orbitals):


    data_dict[str(atom.orbital_info.n_qnums[i])+"_"+str(atom.orbital_info.k_qnums[i])+'_large_div_r']=atom.orbital_info.interpolated_ldivr_components[i]
    data_dict[str(atom.orbital_info.n_qnums[i])+"_"+str(atom.orbital_info.k_qnums[i])+'_nlarge_div_r']=atom.orbital_info.interpolated_nldivr_components[i]
    data_dict[str(atom.orbital_info.n_qnums[i])+"_"+str(atom.orbital_info.k_qnums[i])+'_small_div_r']=atom.orbital_info.interpolated_sdivr_components[i]
    data_dict[str(atom.orbital_info.n_qnums[i])+"_"+str(atom.orbital_info.k_qnums[i])+'_nsmall_div_r']=atom.orbital_info.interpolated_nsdivr_components[i]

  IP_Data=pd.DataFrame.from_dict(data_dict)
  IP_Data.to_csv('interpolated_orbitals.csv')
 
  print('Time to run: '+str(time.time()-stime)+'s.')
