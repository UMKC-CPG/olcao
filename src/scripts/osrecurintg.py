#!/usr/bin/env python3

import argparse as ap
import os
import sys
import osrecurintgana as ana
import osrecurintgnum as num

"""
PROGRAM: osrecurintg.py
PURPOSE: This program implements the Obara-Saika recursive method for computing
         one-electron integrals. The output is Fortran source code designed to
         be used in the OLCAO program suite or a separate testing program.
"""

class ScriptSettings():
    def __init__(self, temp_params):

        # Output format.
        self.production = temp_params[0]

        # Whether or not to vectorize the production output.
        self.vectorize = temp_params[1]

        # Highest orbital (l) angular momentum.
        self.max_lam = temp_params[2]

        # List of integrals to produce.
        self.overlap = temp_params[3]
        self.kinetic = temp_params[4]
        self.nuclear = temp_params[5]
        self.electron = temp_params[6]
        self.momentum = temp_params[7]


def init_environment():
    """Define default values for the integral equations by pulling them from
    the resource control file $OLCAO_DIR/.olcao/osrecurintgrc.py."""

    # Read default variables from the resource control file.
    sys.path.insert(1, os.getenv('OLCAO_DIR') + "/.olcao")
    from osrecurintgrc import parameters_and_defaults

    # Create and return an object containing the script settings.
    return ScriptSettings(parameters_and_defaults())


def parse_command_line():

    # Create the parser tool.
    parser = ap.ArgumentParser(description='Control parameters')

    # Add arguments to the parser.
    parser.add_argument('-p', '--production', action='store_true',
                        help='Request that code be produced for use in a '+
                        'production environment instead of the default '+
                        'behavior of producing code for use in a testing '+
                        'environment.', default='store_false',)

    parser.add_argument('-v', '--vectorize', action='store_true',
                        default=False, help = 'Vectorize the '+
                        'integrals. (Option ignored if not also production.)')

    parser.add_argument('-l', '--max_lam', type=int, help='Highest angular '+
                        'momentum to produce equations for.', default=-1)

    parser.add_argument('-o', '--overlap', action='store_true',
                        default=False, help = 'Include two-center '+
                        'overlap integral.')

    parser.add_argument('-k', '--kinetic', action='store_true',
                        default=False, help = 'Include two-center '+
                        'kinetic energy integral.')

    parser.add_argument('-n', '--nuclear', action='store_true',
                        default=False, help = 'Include three-center '+
                        'nuclear attraction integral.')

    parser.add_argument('-e', '--electron', action='store_true',
                        default=False, help = 'Include three-center '+
                        'overlap integral (with third term being s-type), '+
                        'which is used for electron repulsion.')

    parser.add_argument('-m', '--momentum', action='store_true',
                        default=False, help = 'Include two-center '+
                        'momentum integral.')

    parser.add_argument('-a', '--all', action='store_true',
                        default=False, help = 'Activate o,k,n,e,m '+
                        'integral types.')

    return parser.parse_args()


def reconcile(args, settings):

    # Copy each given argument over the default settings. If an argument was
    #   not given, then retain the default setting for that argument.
    if (args.production == True):
        settings.production = True

    if (args.vectorize == True):
        settings.vectorize = True

    if (args.max_lam != -1):
        settings.max_lam = args.max_lam

    if (args.overlap == True):
        settings.overlap = True

    if (args.kinetic == True):
        settings.kinetic = True

    if (args.nuclear == True):
        settings.nuclear = True

    if (args.electron == True):
        settings.electron = True

    if (args.momentum == True):
        settings.momentum = True

    if (args.all == True):
        settings.overlap = True
        settings.kinetic = True
        settings.nuclear = True
        settings.electron = True
        settings.momentum = True

    # Return settings that incorporate the command line arguments.
    return settings


def create_data_aids(max_lam):

    # Note, for the conversion matrix each row has two groups of three terms
    #   each, which implies that each solution matrix element for the integral
    #   of spherical-harmonic GTOs can be composed of up to three terms from
    #   the matrix of primitive Cartesian GTO integral solutions.

    # Refer to Molecular Electronic Structure Theory (Chapter 9) for the
    #   definitions of spherical harmonic GTOs and primitive Cartesian GTOs.
    # The pc-GTOs follow: G_ijk = x^i * y^j * z^k * exp(-alpha*r^2).
    # The sh-GTOs follow: G_lm = Y_lm * exp(-alpha*r^2) which can be
    #   expressed in terms of pc-GTOs. (The Y_lm are real spherical
    #   harmonics.)

    # The conversion values in the first group of three are coefficients.
    #   The values in the second groups are the indices within the pc matrix
    #   to pull from to create the sh matrix. In the future, this should
    #   probably all be reordered to either align with the form given in
    #   equation 9.1.9 of Molecular Electronic Structure Theory, or it
    #   should be ordered to optimize computational performance.

    # Here is a reference for the sh matrix out to f orbitals as ordered in
    #   the OLCAO method (and encoded in the conversion matrix):
    # sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
    # 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
    # 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy
 
    # Here is a reference for the pc matrix out to f orbitals as ordered in
    #   the OLCAO method (and encoded in the conversion matrix):
    # pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
    # 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
    # 18,xxx; 19,yyy; 20,zzz

    # Add the s conversion, s triad, and l angular momentum total for sh, pc.
    conversion = [[[1, 0, 0], [1, 1, 1]]]  # 1 s
    triads = [[0, 0, 0]] # s 0
    lam_sh_list = [0]
    lam_pc_list = [0]

    # Compute the number of spherical harmonic GTOs.
    num_sh_intg = int((2*max_lam + 2)*(max_lam + 1) / 2)
    
    # Begin accumulating the number of primitive Cartesian GTOs. (Note "=")
    num_pc_intg = int((0 + 2)*(0 + 1) / 2)

    if (max_lam == 0):
        return conversion, triads, lam_sh_list, lam_pc_list, num_sh_intg, \
                num_pc_intg

    # Add the p conversions, p triads, and l angular momentum totals.
    conversion.append([[1, 0, 0], [2, 1, 1]])  # 2 p x
    conversion.append([[1, 0, 0], [3, 1, 1]])  # 3 p y
    conversion.append([[1, 0, 0], [4, 1, 1]])  # 4 p z
    triads.append([1, 0, 0]) # x 1
    triads.append([0, 1, 0]) # y 2
    triads.append([0, 0, 1]) # z 3
    for i in range(3):
        lam_sh_list.append(1)
        lam_pc_list.append(1)

    # Accumulate an additional number of primitive Cartesian GTOs. (Note "+=")
    num_pc_intg += int((1 + 2)*(1 + 1) / 2)

    if (max_lam == 1):
        return conversion, triads, lam_sh_list, lam_pc_list, num_sh_intg, \
                num_pc_intg

    # Add the d conversions, d triads, and l angular momentum sh, pc totals.
    conversion.append([[1, 0, 0], [8, 1, 1]])    # 5 d xy
    conversion.append([[1, 0, 0], [9, 1, 1]])    # 6 d xz
    conversion.append([[1, 0, 0], [10, 1, 1]])   # 7 d yz
    conversion.append([[1, -1, 0], [5, 6, 1]])   # 8 d xx-yy
    conversion.append([[2, -1, -1], [7, 5, 6]])  # 9 d 2zz-xx-yy
    triads.append([2, 0, 0]) # xx 4
    triads.append([0, 2, 0]) # yy 5
    triads.append([0, 0, 2]) # zz 6
    triads.append([1, 1, 0]) # xy 7
    triads.append([1, 0, 1]) # xz 8
    triads.append([0, 1, 1]) # yz 9
    for i in range(5):
        lam_sh_list.append(2)
    for i in range(6):
        lam_pc_list.append(2)

    # Accumulate an additional number of primitive Cartesian GTOs. (Note "+=")
    num_pc_intg += int((2 + 2)*(2 + 1) / 2)

    if (max_lam == 2):
        return conversion, triads, lam_sh_list, lam_pc_list, num_sh_intg, \
                num_pc_intg

    # Add the f conversions, f triads, and l angular momentum totals.
    conversion.append([[1, 0, 0], [11, 1, 1]])      # 10 f xyz
    conversion.append([[1, -1, 0], [13, 15, 1]])    # 11 f xxz-yyz
    conversion.append([[1, -3, 0], [18, 14, 1]])    # 12 f xxx-3yyx
    conversion.append([[3, -1, 0], [12, 19, 1]])    # 13 f 3xxy-yyy
    conversion.append([[2, -3, -3], [20, 13, 15]])  # 14 f 2zzz-3xxz-3yyz
    conversion.append([[4, -1, -1], [16, 18, 14]])  # 15 f 4zzx-xxx-yyx
    conversion.append([[4, -1, -1], [17, 12, 19]])  # 16 f 4zzy-xxy-yyy
    triads.append([1, 1, 1]) # xyz 10
    triads.append([2, 1, 0]) # xxy 11
    triads.append([2, 0, 1]) # xxz 12
    triads.append([1, 2, 0]) # xyy 13
    triads.append([0, 2, 1]) # yyz 14
    triads.append([1, 0, 2]) # xzz 15
    triads.append([0, 1, 2]) # zzy 16
    triads.append([3, 0, 0]) # xxx 17
    triads.append([0, 3, 0]) # yyy 18
    triads.append([0, 0, 3]) # zzz 19
    for i in range(7):
        lam_sh_list.append(3)
    for i in range(10):
        lam_pc_list.append(3)

    # Accumulate an additional number of primitive Cartesian GTOs. (Note "+=")
    num_pc_intg += int((3 + 2)*(3 + 1) / 2)

    if (max_lam == 3):
        return conversion, triads, lam_sh_list, lam_pc_list, num_sh_intg, \
                num_pc_intg

    # Add the sh conversions, sh triads, and l angular momentum totals.
    triads.append([0, 1, 3]) # yzzz
    triads.append([0, 3, 1]) # yyyz
    triads.append([1, 0, 3]) # xzzz
    triads.append([1, 3, 0]) # xyyy
    triads.append([3, 0, 1]) # xxxz
    triads.append([3, 1, 0]) # xxxy
    triads.append([4, 0, 0]) # xxxx
    triads.append([0, 4, 0]) # yyyy
    triads.append([0, 0, 4]) # zzzz
    triads.append([1, 1, 2]) # xyzz
    triads.append([1, 2, 1]) # xyyz
    triads.append([2, 1, 1]) # xxyz
    for i in range(9): # CHECK
        lam_sh_list.append(4)
    for i in range(12): # CHECK
        lam_pc_list.append(4)

    # Accumulate an additional number of primitive Cartesian GTOs. (Note "+=")
    num_pc_intg += int((4 + 2)*(4 + 1) / 2)

    if (max_lam == 4):
        return conversion, triads, lam_sh_list, lam_pc_list, num_sh_intg, \
                num_pc_intg


def print_head(settings, f, num_sh_intg, num_pc_intg, max_lam):

    if (settings.production):
        head = """module O_GaussianIntegralsVec"""
    else:
        head = """program GaussianIntegrals"""

    head += """

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none"""

    if (settings.production):
        head += """

   ! Define the access level.
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

   ! These subroutines are based primarily on the work of Obara-Saika as
   !   published in: Obara S, Saika A., "Efficient recursive computation of
   !   molecular integrals over Cartesian Gaussian functions", The Journal of
   !   Chemical Physics, (1986), Volume 84(7), pages: 3963-3974. Additional
   !   guidance can be found in: Helgaker T, Jørgensen P, Olsen J.,
   !   "Molecular electronic-structure theory. Chichester ; New York: Wiley;
   !   2000. 908 pages; (Primarily chapter 9).
   
"""

    if (not settings.production):
        head += """

   ! Define program variables.
   real (kind=double), allocatable, dimension (:,:,:,:) :: pc
   real (kind=double), allocatable, dimension (:,:,:,:) :: sh
   real (kind=double), dimension (3) :: alphas
   real (kind=double), dimension (3,3) :: pos ! x,y,z of A,B,C
   real (kind=double), dimension (3,2) :: temp_alphas
   real (kind=double), dimension (3,3,2) :: temp_pos ! x,y,z of A,B,C
   real (kind=double), dimension (3) :: temp_alphas_step
   real (kind=double), dimension (3,3) :: temp_pos_step ! x,y,z of A,B,C
   real (kind=double) :: cell_size, step_size ! Parameters for the numerical
         ! integration process. Read carefully! The cell size is given as a
         ! positive number that defines the maximum extent of the space to
         ! integrate along any given axis. Hence the actual cell is double
         ! that number. I.e., integrate along x, y, z axes from -cell_size
         ! to +cell_size. The step size is just what it appears to be, but
         ! note that it is for numerical integration and is a very different
         ! concept and is unrelated to the "num_steps" defined below.
   integer :: num_segments, num_steps
   integer :: h, i, p, q

   ! Allocate space to hold the appropriately sized pc and sh matrices. The
   !   last index in both pc and sh is a 2 to hold analytical solutions
   !   (value = 1) and numerical solutions (value = 2). The second to last
   !   index is the max of three or (2*max_lam+1) to allow for the momentum
   !   matrix and multi-level nuclear solutions (but all elements are not
   !   always used). The pc matrix is used multiple times to produce the three
   !   momentum matrix sets (x,y,z).
   allocate (pc(""" \
    + f"{num_pc_intg},{num_pc_intg},{max(2*max_lam+1,3)},2))" + """
   allocate (sh(""" + f"{num_sh_intg},{num_sh_intg},3,2))" + """

   ! Open the control file.
   open (10, file="intgcontrol", status='old')

   ! The first bits of information in the control file will control the
   !   numerical integration. The
   read (10,*) cell_size, step_size

   ! Read the number of segments.
   read (10,*) num_segments

   ! For each segment, read in the number of steps, allocate space to hold
   !   relevant information for each step, compute the step descriptions,
   !   compute the integral matrices, store the results, deallocate for the
   !   next segment.
   do h = 1, num_segments

      ! Read in the number of steps for this segment.
      read (10,*) num_steps
      if (num_steps == 1) then
   
         ! Read the parameters for the current segment.
         read (10,*) temp_alphas(:,1), temp_pos(:,:,1)

         ! Fill in dummy values for the temp step data. (It will not be used.)
         temp_alphas_step(:) = 0.0d0
         temp_pos_step(:,:) = 0.0d0
      else
   
         ! Read the beginning and ending parameters for the current segment.
         read (10,*) temp_alphas(:,1), temp_pos(:,:,1) ! First step
         read (10,*) temp_alphas(:,2), temp_pos(:,:,2) ! Last step

         ! Compute the size of the step for each parameter.
         temp_alphas_step(:) = (temp_alphas(:,2) - temp_alphas(:,1))/num_steps
         temp_pos_step(:,:) = (temp_pos(:,:,2) - temp_pos(:,:,1))/num_steps
      endif

      ! Start iterating over the steps.
      do i = 1, num_steps
         alphas(:) = temp_alphas(:,1) + (i-1) * temp_alphas_step(:)
         pos(:,:) = temp_pos(:,:,1) + (i-1) * temp_pos_step(:,:)
"""
        if (settings.overlap):
            head += """

         ! Compute the pc and sh integral results for the current parameters
         !   using the analytic formulas.
         call overlap2CIntgAna(alphas(1),alphas(2),pos(:,1),pos(:,2),&
               & pc(:,:,1,1),sh(:,:,1,1))

         ! Compute the pc and sh integral results for the current parameters
         !   using numerical integration.
         call overlap2CIntgNum(alphas(1),alphas(2),pos(:,1),pos(:,2),&
               & pc(:,:,1,2),sh(:,:,1,2),cell_size,step_size)

         ! Print the pc and sh integral result differences.
         call print_pc_sh(h,i,1,alphas,pos,pc(:,:,1,:),sh(:,:,1,:),&
               & "overlap.dat")
"""

        if(settings.kinetic):
            head += """

         ! Compute the pc and sh integral results for the current parameters
         !   using analytical formulas.
         call kinetic2CIntgAna(alphas(1),alphas(2),pos(:,1),pos(:,2),&
               & pc(:,:,1,1),sh(:,:,1,1))

         ! Compute the pc and sh integral results for the current parameters
         !   using numerical integration.
         call kinetic2CIntgNum(alphas(1),alphas(2),pos(:,1),pos(:,2),&
               & pc(:,:,1,2),sh(:,:,1,2),cell_size,step_size)

         ! Print the pc and sh integral result differences.
         call print_pc_sh(h,i,2,alphas,pos,pc(:,:,1,:),sh(:,:,1,:),&
               & "kinetic.dat")
"""

        if(settings.nuclear):
            head += """

         ! Compute the pc and sh integral results for the current parameters
         !   using analytical formulas.
         call nuclear3CIntgAna(alphas(1),alphas(2),alphas(3),pos(:,1),&
               & pos(:,2),pos(:,3),pc(:,:,:,1),sh(:,:,1,1))

         ! Compute the pc and sh integral results for the current parameters
         !   using numerical integration.
         call nuclear3CIntgNum(alphas(1),alphas(2),alphas(3),pos(:,1),&
               & pos(:,2),pos(:,3),pc(:,:,1,2),sh(:,:,1,2),cell_size,&
               & step_size)

         ! Print the pc and sh integral result differences.
         call print_pc_sh(h,i,3,alphas,pos,pc(:,:,1,:),sh(:,:,1,:),&
               & "nuclear.dat")
"""

        if(settings.electron):
            head += """

         ! Compute the pc and sh integral results for the current parameters
         !   using analytical formulas.
         call electron3CIntgAna(alphas(1),alphas(2),alphas(3),pos(:,1),&
               & pos(:,2),pos(:,3),pc(:,:,1,1),sh(:,:,1,1))

         ! Compute the pc and sh integral results for the current parameters
         !   using numerical integration.
         call electron3CIntgNum(alphas(1),alphas(2),alphas(3),pos(:,1),&
               & pos(:,2),pos(:,3),pc(:,:,1,2),sh(:,:,1,2),cell_size,&
               & step_size)

         ! Print the pc and sh integral result differences.
         call print_pc_sh(h,i,4,alphas,pos,pc(:,:,1,:),sh(:,:,1,:),&
               & "electrn.dat")
"""

        if(settings.momentum):
            head += """

         ! Compute the pc and sh integral results for the current parameters.
         call momentum2CIntgAna(alphas(1),alphas(2),pos(:,1),pos(:,2),&
               & """
            head += f"pc(:,:,:,1),sh(:,:,:,1))"
            head += """

         ! Compute the pc and sh integral results for the current parameters
         !   using numerical integration.
         call momentum2CIntgNum(alphas(1),alphas(2),pos(:,1),pos(:,2),&
               & """
            head += f"pc(:,:,:,2),sh(:,:,:,2),cell_size,"
            head += """step_size)

         ! Print the pc and sh integral result differences.
"""

            for xyz in range(3):
                head += "         "
                head += f"call print_pc_sh(h,i,{xyz+5},alphas,pos,"
                head += f"pc(:,:,{xyz+1},:),sh(:,:,{xyz+1},:),&\n"
                if (xyz == 0):
                    head += '               & "momentx.dat")\n'
                elif (xyz == 1):
                    head += '               & "momenty.dat")\n'
                elif (xyz == 2):
                    head += '               & "momentz.dat")\n'
        head += """
      enddo
   enddo

   ! Deallocate space.
   deallocate(pc)
   deallocate(sh)

"""

    # Print the one subroutine shared by all other subroutines.
        head += """
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module/program subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

   subroutine print_pc_sh(h,i,unitPlus,alphas,pos,pc,sh,filename)

      ! Use necessary modules
      use O_Kinds

      ! Make sure no funny variables are created.
      implicit none

      ! Define dummy variables.
      integer :: h, i, unitPlus
      real(kind=double), dimension(3) :: alphas
      real(kind=double), dimension(3,3) :: pos ! xyz (1st idx) of ABC (2nd)
      real(kind=double), dimension("""
        head += f"{num_pc_intg},{num_pc_intg},2) :: pc"
        head += """ ! Last idx 1 = ana; 2 = num
      real(kind=double), dimension("""
        head += f"{num_sh_intg},{num_sh_intg},2) :: sh"
        head += """ ! Last idx 1 = ana; 2 = num
      character*11 :: filename

      ! Define local variables.
      logical :: io_opened
      integer :: j,k

      ! Open the output file if it isn't already open.
      inquire(299+unitPlus, OPENED = io_opened)
      if (.not. io_opened) then
         open (299+unitPlus, file=filename, status="unknown")
      endif

      write (299+unitPlus,fmt="(2i5)",advance="NO") h, i
      write (299+unitPlus,fmt="(3e16.8)",advance="NO") alphas(:)
      write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,1)
      write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,2)
      write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,3)

      ! Print the result difference for the given pc matrix.
      do j = 1, """ + f"{num_pc_intg}" + """
         do k = 1, """ + f"{num_pc_intg}" + """
            write (299+unitPlus,fmt="(e16.8)",advance="NO") pc(k,j,1)
            write (299+unitPlus,fmt="(e16.8)",advance="NO") pc(k,j,2)
            write (299+unitPlus,fmt="(e16.8)",advance="NO") &
                  & pc(k,j,1) - pc(k,j,2)
         enddo
      enddo

      ! Print the result difference for the given sh matrix.
      do j = 1, """ + f"{num_sh_intg}" + """
         do k = 1, """ + f"{num_sh_intg}" + """
            write (299+unitPlus,fmt="(e16.8)",advance="NO") &
                  & sh(k,j,1) - sh(k,j,2)
         enddo
      enddo

      ! Write an endline to prepare for the next step.
      write (299+unitPlus, *) ""
   end subroutine
"""

    f.write(head)


def print_foot(production, f):
    if (production):
        foot = """
end module O_GaussianIntegralsVec
"""
    else:
        foot = """
end program GaussianIntegrals
"""

    f.write(foot)


def main():

    # Set default values from the user resource control file.
    settings = init_environment()

    # Parse the command line.
    args = parse_command_line()

    # Reconcile the default settings and the given arguments
    settings = reconcile(args, settings)

    # Create the data aids for the requested set of angular momenta. This
    #   includes a matrix to aid the conversion from the pc matrix to the sh
    #   matrix, the triads defining the Cartesian indices of each orbital,
    #   the total l angular momentum of each triad, and the number of terms
    #   in each of the matrices to produce (num_sh_intg for the sh-matrix
    #   and num_pc_intg for the pc-matrix).
    conversion, triads, lam_sh_list, lam_pc_list, num_sh_intg, num_pc_intg \
            = create_data_aids(settings.max_lam)

    # The basic structure of the program is to allow for two different output
    #   types: production and test.
    # In the production environment we will create a module that contains one
    #   subroutine for each type of integral. Each subroutine will produce a
    #   sh-matrix result for the given set of alphas, site coordinates, and
    #   the specific choice of the l1l2switch.
    # In the test environment we will create two subroutines for each type of
    #   integral. The first subroutine will produce a sh-matrix and pc-matrix
    #   result for the given set of alphas and site coordinates. The
    #   additional difference from the production environment is that the
    #   test subroutines will always return full results for each matrix.
    #   I.e., like a maximal l1l2switch call. The second subroutine in the
    #   test environment will compute the values for the full sh and pc
    #   matrices using a numerical algorithm. The final difference between
    #   the test and production environment is that the test environment will
    #   include code to read control parameters for the integrals that will
    #   actually call the integral subroutines and print comparisons of their
    #   results.

    # Open the output file for the program.
    if (settings.production):
        if (settings.vectorize):
            f = open("gaussIntegrals.vec.f90", "w")
        else:
            f = open("gaussIntegrals.f90", "w")
    else:
        f = open("gaussIntgTest.f90", "w")

    # Print the overall head.
    print_head(settings, f, num_sh_intg, num_pc_intg, settings.max_lam)

    # Proceed with computing and printing the results.

    # Manage the two-center overlap integrals.
    if (settings.overlap):
        matrix = ana.overlap(triads, False, False, settings.vectorize)
        if (settings.production):
            if (settings.vectorize):
                ana.print_production_overlap_vec(conversion, triads, matrix, 
                                                 lam_sh_list, lam_pc_list, f)
            else:
                ana.print_production_overlap(conversion, triads, matrix,
                                             lam_sh_list, lam_pc_list,f)
        else:
            ana.print_test_overlap_ana(conversion, triads, matrix, f)
            num.print_test_overlap_num(conversion, triads, f)

    # Manage the kinetic energy integrals.
    if (settings.kinetic):
        matrix_ol = ana.overlap(triads, True, False, settings.vectorize)
        matrix_ke = ana.kinetic(triads, matrix_ol, settings.vectorize)
        if (settings.production):
            if (settings.vectorize):
                ana.print_production_kinetic_vec(conversion, triads, matrix_ke,
                                                 matrix_ol, lam_sh_list,
                                                 lam_pc_list, f)
            else:
                ana.print_production_kinetic(conversion, triads, matrix_ke,
                                             matrix_ol, lam_sh_list,
                                             lam_pc_list, f)
        else:
            ana.print_test_kinetic_ana(conversion, triads, matrix_ke,
                                       matrix_ol, f)
            num.print_test_kinetic_num(conversion, triads, f)

    # Manage the three-center overlap (e-e repulsion) integrals.
    if (settings.electron):
        # The recursive matrix terms for the three center integral are all
        #   exactly the same as for the two center. The only differences
        #   will be evident in the definition of some base variables in the
        #   printed Fortran subroutine. (I.e., zeta, and G instead of P.)
        matrix = ana.overlap(triads, False, True, settings.vectorize)
        if (settings.production):
            if (settings.vectorize):
                ana.print_production_electron_vec(conversion, triads, matrix,
                                                  lam_sh_list, lam_pc_list, f)
            else:
                ana.print_production_electron(conversion, triads, matrix,
                                              lam_sh_list, lam_pc_list,f)
        else:
            ana.print_test_electron_ana(conversion, triads, matrix, f)
            num.print_test_electron_num(conversion, triads, f)

    # Manage the nuclear attraction potential integrals.
    if (settings.nuclear):
        num_triads = len(triads)
        matrix = [[[""]*(2*settings.max_lam+1) for i in range(num_triads)] \
                   for j in range(num_triads)]
        for m in range(2*settings.max_lam + 1):
            matrix_temp = ana.nuclear(triads, m, settings.vectorize)
            for i in range(num_triads):
                for j in range(num_triads):
                    matrix[i][j][m] = matrix_temp[i][j]
        if (settings.production):
            if (settings.vectorize):
                ana.print_boys_vec(f)
                ana.print_production_nuclear_vec(conversion, triads, matrix,
                                                 settings.max_lam,
                                                 lam_sh_list, lam_pc_list, f)
            else:
                ana.print_boys(f)
                ana.print_production_nuclear(conversion, triads, matrix,
                                             settings.max_lam, lam_sh_list,
                                             lam_pc_list, f)
        else:
            ana.print_test_nuclear_ana(conversion, triads, matrix,
                                       settings.max_lam, f)
            ana.print_boys(f)
            num.print_test_nuclear_num(conversion, triads, f)

    # Manage the momentum matrix integrals.
    if (settings.momentum):
        matrix_ol = ana.overlap(triads, True, False, settings.vectorize)
        matrix_mm = ana.momentum(triads, settings.vectorize)
        if (settings.production):
            if (settings.vectorize):
                ana.print_production_momentum_vec(conversion, triads,
                                                  matrix_mm, matrix_ol,
                                                  lam_sh_list, lam_pc_list, f)
            else:
                ana.print_production_momentum(conversion, triads, matrix_mm,
                                              matrix_ol, lam_sh_list,
                                              lam_pc_list, f)
        else:
            ana.print_test_momentum_ana(conversion, triads, matrix_mm,
                                        matrix_ol, f)
            num.print_test_momentum_num(conversion, triads, f)

    # Print the overall foot.
    print_foot(settings.production, f)

    # Close the program output file.
    f.close()


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or a request
    #   to import information from external modules. Only now do we actually
    #   start running the program.
    main()
