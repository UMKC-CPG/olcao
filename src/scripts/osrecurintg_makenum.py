#!/usr/bin/env python3

"""
PROGRAM: osrecurintg_makenum.py
PURPOSE: This program is used to produce the equations that will be used to
produce numerical solutions to multicenter gaussian type orbital integrals.
I.e., this program prints out chunks of Fortran code that can be directly
copied and pasted into other python programs that are used to create actual
Fortran programs for computing numerical solutions. (The reason we are
separating the creation of the equations from the creation of the programs
is that the creation of the equations can be somewhat time consuming to run
and once they have been determined, there is really no good reason to
recreate them each time we do a test (or production) run of the Obara-Saika
programs.)

Although the output equations may be directly copied, that does not mean that
it is the best way to do things. The equations that are produced might be
separated into distinct subcomponents that can be organized into nested loops
to increase the computational efficiency. Depending on the integration mesh
that is used, the numerical integral computation may be time consuming.
Therefore, it is probably worthwhile to put some effort into reorganizing
the equation subcomponents. That being said, it is absolutely essential that
every effort is made to ensure that no errors are introduced when performing
the reorganization. If there is a test that can be performed to ensure that
no errors are introduced, then please try to do that test.
"""

import sympy as sp
import argparse as ap
import os
import sys
import re
import osrecurintglib as lib


class ScriptSettings():
    def __init__(self, temp_params):

        # List of integrals to produce.
        self.overlap = temp_params[0]
        self.kinetic = temp_params[1]
        self.nuclear = temp_params[2]
        self.electron = temp_params[3]
        self.momentum = temp_params[4]
        self.massvel = temp_params[5]
        self.dipole = temp_params[6]
        self.dkinetic = temp_params[7]
        self.dnuclear = temp_params[8]
        self.delectron = temp_params[9]
        # self.e_field = temp_params[5]
        # self.e_grad = temp_params[6]
        # self.ang_mom = temp_params[7]
        # self.spin_orb = temp_params[8]


def init_environment():
    """Define default values for the integral equations by pulling them from
    the resource control file $OLCAO_DIR/.olcao/osrecurintgrc_makenum.py."""

    # Read default variables from the resource control file.
    sys.path.insert(1, os.getenv('OLCAO_DIR') + "/.olcao")
    from osrecurintg_makenumrc import parameters_and_defaults

    # Create and return an object containing the script settings.
    return ScriptSettings(parameters_and_defaults())


def parse_command_line():

    # Create the parser tool.
    parser = ap.ArgumentParser(description='Control parameters')

    # Add arguments to the parser.
    parser.add_argument('-o', '--overlap', action='store_true',
                        default='store_false', help = 'Include two-center '+
                        'overlap integral.')

    parser.add_argument('-k', '--kinetic', action='store_true',
                        default='store_false', help = 'Include two-center '+
                        'kinetic energy integral.')

    parser.add_argument('-n', '--nuclear', action='store_true',
                        default='store_false', help = 'Include three-center '+
                        'nuclear attraction integral.')

    parser.add_argument('-e', '--electron', action='store_true',
                        default='store_false', help = 'Include three-center '+
                        'overlap integral (with third term being s-type), '+
                        'which is used for electron repulsion.')

    parser.add_argument('-m', '--momentum', action='store_true',
                        default='store_false', help = 'Include two-center '+
                        'momentum integral.')

    parser.add_argument('-mv', '--massvel', action='store_true',
                        default='store_false', help = 'Include mass '+
                        'velocity integral.')

    parser.add_argument('-d', '--dipole', action='store_true',
                        default='store_false', help = 'Include the dipole '+
                        'dipole integral.')

    parser.add_argument('-dk', '--dkinetic', action='store_true',
                        default='store_false', help = 'Include two-center '+
                        'derivative of the kinetic energy integral.')

    parser.add_argument('-dn', '--dnuclear', action='store_true',
                        default='store_false', help = 'Include three-center '+
                        'derivative of the nuclear attraction integral.')

    parser.add_argument('-de', '--delectron', action='store_true',
                        default='store_false', help = 'Include three-center '+
                        'overlap integral (3rd term being s-type), which is '+
                        'used for the derivative of the electron repulsion.')

    return parser.parse_args()


def reconcile(args, settings):

    # Copy each given argument over the default settings. If an argument was
    #   not given, then retain the default setting for that argument.
    
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

    if (args.massvel == True):
        settings.massvel = True

    if (args.dipole == True):
        settings.dipole = True

    if (args.dkinetic == True):
        settings.dkinetic = True

    if (args.dnuclear == True):
        settings.dnuclear = True

    if (args.delectron == True):
        settings.delectron = True

    # Return settings that incorporate the command line arguments.
    return settings

# This integral is seperable. We only need to solve it for one dimension and
#   then we will use string substitutions to extend it to three dimensions.
def overlap():
    Px = sp.symbols('Px')
    lx1 = sp.symbols('lx1')
    lx2 = sp.symbols('lx2')
    zeta1, zeta2 = sp.symbols('zeta1 zeta2')
    Rx1 = sp.symbols('Rx1')
    Rx2 = sp.symbols('Rx2')

    orbital_1 = (Px-Rx1)**lx1 * sp.exp(-zeta1*(Px-Rx1)**2)
    orbital_2 = (Px-Rx2)**lx2 * sp.exp(-zeta2*(Px-Rx2)**2)

    return sp.printing.fcode(orbital_1 * orbital_2)


def kinetic_energy():
    # Px = sp.symbols('Px')
    # lx1 = sp.symbols('lx1')
    # lx2 = sp.symbols('lx2')
    # zeta1, zeta2 = sp.symbols('zeta1 zeta2')
    # Rx1 = sp.symbols('Rx1')
    # Rx2 = sp.symbols('Rx2')

    # orbital_1 = (Px-Rx1)**lx1 * sp.exp(-zeta1*(Px-Rx1)**2)
    # orbital_2 = (Px-Rx2)**lx2 * sp.exp(-zeta2*(Px-Rx2)**2)

    # # Note, the factor of -1/2 is accounted for in the osrecurint.py program.
    # return sp.printing.fcode(orbital_1 * sp.diff(sp.diff(orbital_2, Px), Px))

    # The commented code below will produce the full expression. However,
    #   the uncommented code above will produce the expression for one
    #   dimension. That 1D expression is what we later manipulate by hand
    #   to create the separable solution.

    # Px, Py, Pz = sp.symbols('Px Py Pz')
    # lx1, ly1, lz1 = sp.symbols('lx1 ly1 lz1')
    # lx2, ly2, lz2 = sp.symbols('lx2 ly2 lz2')
    # zeta1, zeta2, zeta3 = sp.symbols('zeta1 zeta2 zeta3')
    # Rx1, Ry1, Rz1 = sp.symbols('Rx1 Ry1 Rz1')
    # Rx2, Ry2, Rz2 = sp.symbols('Rx2 Ry2 Rz2')
    # Rx3, Ry3, Rz3 = sp.symbols('Rx3 Ry3 Rz3')

    # orbital_1 = (Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * \
    #         sp.exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))
    # orbital_2 = (Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 * \
    #         sp.exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))

    # return sp.printing.fcode(
    #         orbital_1 * sp.diff(sp.diff(orbital_2, Px), Px)
    #         + orbital_1 * sp.diff(sp.diff(orbital_2, Py), Py)
    #         + orbital_1 * sp.diff(sp.diff(orbital_2, Pz), Pz))

    # (s|s) Case:
    Px, Py, Pz = sp.symbols('Px Py Pz')
    lx1, ly1, lz1 = sp.symbols('lx1 ly1 lz1')
    lx2, ly2, lz2 = sp.symbols('lx2 ly2 lz2')
    zeta1, zeta2, zeta3 = sp.symbols('zeta1 zeta2 zeta3')
    Rx1, Ry1, Rz1 = sp.symbols('Rx1 Ry1 Rz1')
    Rx2, Ry2, Rz2 = sp.symbols('Rx2 Ry2 Rz2')
    Rx3, Ry3, Rz3 = sp.symbols('Rx3 Ry3 Rz3')

    orbital_1 = sp.exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))
    orbital_2 = sp.exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))

    return sp.printing.fcode(
            orbital_1 * sp.diff(sp.diff(orbital_2, Px), Px)
            + orbital_1 * sp.diff(sp.diff(orbital_2, Py), Py)
            + orbital_1 * sp.diff(sp.diff(orbital_2, Pz), Pz))


# This integral is /not/ seperable and so we must produce the full 3D
#   solution.
def nuclear():
    Px, Py, Pz = sp.symbols('Px Py Pz')
    lx1, ly1, lz1 = sp.symbols('lx1 ly1 lz1')
    lx2, ly2, lz2 = sp.symbols('lx2 ly2 lz2')
    zeta1, zeta2, zeta3 = sp.symbols('zeta1 zeta2 zeta3')
    Rx1, Ry1, Rz1 = sp.symbols('Rx1 Ry1 Rz1')
    Rx2, Ry2, Rz2 = sp.symbols('Rx2 Ry2 Rz2')
    Rx3, Ry3, Rz3 = sp.symbols('Rx3 Ry3 Rz3')

    orbital_1 = (Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * \
            sp.exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))
    orbital_2 = (Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 * \
            sp.exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))
    orbital_3 = sp.exp(-zeta3*( (Px-Rx3)**2 + (Py-Ry3)**2 + (Pz-Rz3)**2 ))

    inv_r3 = 1 / sp.sqrt((Px-Rx3)**2 + (Py-Ry3)**2 + (Pz-Rz3)**2)

    return sp.printing.fcode(orbital_1 * (orbital_3 * inv_r3) * orbital_2)


def electron():
    Px = sp.symbols('Px')
    lx1 = sp.symbols('lx1')
    lx2 = sp.symbols('lx2')
    zeta1, zeta2, zeta3 = sp.symbols('zeta1 zeta2 zeta3')
    Rx1 = sp.symbols('Rx1')
    Rx2 = sp.symbols('Rx2')
    Rx3 = sp.symbols('Rx3')

    orbital_1 = (Px-Rx1)**lx1 * sp.exp(-zeta1*(Px-Rx1)**2)
    orbital_2 = (Px-Rx2)**lx2 * sp.exp(-zeta2*(Px-Rx2)**2)
    orbital_3 = sp.exp(-zeta3*(Px-Rx3)**2) # Always s-type so lx3=0.

    return sp.printing.fcode(orbital_1 * orbital_3 * orbital_2)


def momentum():
    Px = sp.symbols('Px')
    lx1 = sp.symbols('lx1')
    lx2 = sp.symbols('lx2')
    zeta1, zeta2 = sp.symbols('zeta1 zeta2')
    Rx1 = sp.symbols('Rx1')
    Rx2 = sp.symbols('Rx2')

    orbital_1 = (Px-Rx1)**lx1 * sp.exp(-zeta1*(Px-Rx1)**2)
    orbital_2 = (Px-Rx2)**lx2 * sp.exp(-zeta2*(Px-Rx2)**2)

    return sp.printing.fcode(orbital_1 * sp.diff(orbital_2, Px), Px)


def massvel():
    Px, Py = sp.symbols('Px Py')
    lx1, ly1 = sp.symbols('lx1 ly1')
    lx2, ly2 = sp.symbols('lx2 ly2')
    zeta1, zeta2 = sp.symbols('zeta1 zeta2')
    Rx1, Ry1 = sp.symbols('Rx1 Ry1')
    Rx2, Ry2 = sp.symbols('Rx2 Ry2')

    orbital_1a = (Px-Rx1)**lx1 * sp.exp(-zeta1*(Px-Rx1)**2)
    orbital_2a = (Px-Rx2)**lx2 * sp.exp(-zeta2*(Px-Rx2)**2)

    orbital_1b = (Px-Rx1)**lx1 * (Py-Ry1)**ly1 * \
            sp.exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2) )
    orbital_2b = (Px-Rx2)**lx2 * (Py-Ry2)**ly2 * \
            sp.exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2) )

    # Note, any coefficients are accounted for in the osrecurintg.py program.
    return (sp.printing.fcode(orbital_1a *
            sp.diff(sp.diff(sp.diff(sp.diff(orbital_2a, Px), Px), Px), Px)),
            sp.printing.fcode(orbital_1b *
            sp.diff(sp.diff(sp.diff(sp.diff(orbital_2b, Py), Py), Px), Px)))

    # # The commented code below will produce the full expression. However,
    # #   the uncommented code above will produce the expression for one
    # #   dimension. That 1D expression is what we later manipulate by hand
    # #   to create the separable solution for full 3D problem.

    # Px, Py, Pz = sp.symbols('Px Py Pz')
    # lx1, ly1, lz1 = sp.symbols('lx1 ly1 lz1')
    # lx2, ly2, lz2 = sp.symbols('lx2 ly2 lz2')
    # zeta1, zeta2, zeta3 = sp.symbols('zeta1 zeta2 zeta3')
    # Rx1, Ry1, Rz1 = sp.symbols('Rx1 Ry1 Rz1')
    # Rx2, Ry2, Rz2 = sp.symbols('Rx2 Ry2 Rz2')
    # Rx3, Ry3, Rz3 = sp.symbols('Rx3 Ry3 Rz3')

    # orbital_1 = (Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * \
    #         sp.exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))
    # orbital_2 = (Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 * \
    #         sp.exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))

    # # # Specific for the l=0 case:
    # # Px, Py, Pz = sp.symbols('Px Py Pz')
    # # lx1, ly1, lz1 = sp.symbols('lx1 ly1 lz1')
    # # lx2, ly2, lz2 = sp.symbols('lx2 ly2 lz2')
    # # zeta1, zeta2, zeta3 = sp.symbols('zeta1 zeta2 zeta3')
    # # Rx1, Ry1, Rz1 = sp.symbols('Rx1 Ry1 Rz1')
    # # Rx2, Ry2, Rz2 = sp.symbols('Rx2 Ry2 Rz2')
    # # Rx3, Ry3, Rz3 = sp.symbols('Rx3 Ry3 Rz3')

    # # orbital_1 = sp.exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))
    # # orbital_2 = sp.exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))

    # # return sp.printing.fcode(
    # #         orbital_1 * sp.diff(sp.diff(sp.diff(
    # #                             sp.diff(orbital_2, Px), Px), Px), Px)
    # #         + orbital_1 * sp.diff(sp.diff(sp.diff(
    # #                               sp.diff(orbital_2, Py), Py), Py), Py)
    # #         + orbital_1 * sp.diff(sp.diff(sp.diff(
    # #                               sp.diff(orbital_2, Pz), Pz), Pz), Pz)
    # #         + 2 * orbital_1 * sp.diff(sp.diff(sp.diff(
    # #                                  sp.diff(orbital_2, Py), Py), Px), Px)
    # #         + 2 * orbital_1 * sp.diff(sp.diff(sp.diff(
    # #                                  sp.diff(orbital_2, Pz), Pz), Px), Px)
    # #         + 2 * orbital_1 * sp.diff(sp.diff(sp.diff(
    # #                                  sp.diff(orbital_2, Pz), Pz), Py), Py))


def dipole():
        Px, Py, Pz = sp.symbols('Px Py Pz')
        lx1, ly1, lz1 = sp.symbols('lx1 ly1 lz1')
        lx2, ly2, lz2 = sp.symbols('lx2 ly2 lz2')
        mu_x, mu_y, mu_z = sp.symbols('mu_x mu_y mu_z')
        zeta1, zeta2 = sp.symbols('zeta1 zeta2')
        Rx1, Ry1, Rz1 = sp.symbols('Rx1 Ry1 Rz1')
        Rx2, Ry2, Rz2 = sp.symbols('Rx2 Ry2 Rz2')
        Rx3, Ry3, Rz3 = sp.symbols('Rx3 Ry3 Rz3')
        x, y, z = sp.symbols('x y z')

        orbital_1 = (Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * \
                     sp.exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))
        orbital_2 = (Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 * \
                     sp.exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))
        mu = (Px-Rx3)**mu_x * (Py-Ry3)**mu_y * (Pz-Rz3)**mu_z

        return sp.printing.fcode(orbital_1 * mu * orbital_2)


def derivative_kinetic_energy():
    Px, Py = sp.symbols('Px Py')
    lx1, ly1 = sp.symbols('lx1 ly1')
    lx2, ly2 = sp.symbols('lx2 ly2')
    zeta1, zeta2 = sp.symbols('zeta1 zeta2')
    Rx1, Ry1 = sp.symbols('Rx1 Ry1')
    Rx2, Ry2 = sp.symbols('Rx2 Ry2')

    orbital_1a = (Px-Rx1)**lx1 * sp.exp(-zeta1*(Px-Rx1)**2)
    orbital_2a = (Px-Rx2)**lx2 * sp.exp(-zeta2*(Px-Rx2)**2)

    orbital_1b = (Px-Rx1)**lx1 * (Py-Ry1)**ly1 * \
            sp.exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2) )
    orbital_2b = (Px-Rx2)**lx2 * (Py-Ry2)**ly2 * \
            sp.exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2) )

    # Note, any coefficients are accounted for in the osrecurintg.py program.
    return (sp.printing.fcode(orbital_1a *
            sp.diff(sp.diff(sp.diff(orbital_2a, Px), Px), Px)),
            sp.printing.fcode(orbital_1b *
            sp.diff(sp.diff(sp.diff(orbital_2b, Py), Py), Px)))


def derivative_nuclear():
    Px, Py, Pz = sp.symbols('Px Py Pz')
    lx1, ly1, lz1 = sp.symbols('lx1 ly1 lz1')
    lx2, ly2, lz2 = sp.symbols('lx2 ly2 lz2')
    zeta1, zeta2, zeta3 = sp.symbols('zeta1 zeta2 zeta3')
    Rx1, Ry1, Rz1 = sp.symbols('Rx1 Ry1 Rz1')
    Rx2, Ry2, Rz2 = sp.symbols('Rx2 Ry2 Rz2')
    Rx3, Ry3, Rz3 = sp.symbols('Rx3 Ry3 Rz3')

    orbital_1 = (Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * \
            sp.exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))
    orbital_2 = (Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 * \
            sp.exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))
    orbital_3 = sp.exp(-zeta3*( (Px-Rx3)**2 + (Py-Ry3)**2 + (Pz-Rz3)**2 ))

    inv_r3 = 1 / sp.sqrt((Px-Rx3)**2 + (Py-Ry3)**2 + (Pz-Rz3)**2)

    return (sp.printing.fcode(orbital_1 * (
            sp.diff(((orbital_3 * inv_r3) * orbital_2), Px))),
            sp.printing.fcode(orbital_1 * (
            sp.diff(((orbital_3 * inv_r3) * orbital_2), Py))),
            sp.printing.fcode(orbital_1 * (
            sp.diff(((orbital_3 * inv_r3) * orbital_2), Pz))))


def derivative_electron():
    Px = sp.symbols('Px')
    lx1 = sp.symbols('lx1')
    lx2 = sp.symbols('lx2')
    zeta1, zeta2, zeta3 = sp.symbols('zeta1 zeta2 zeta3')
    Rx1 = sp.symbols('Rx1')
    Rx2 = sp.symbols('Rx2')
    Rx3 = sp.symbols('Rx3')

    orbital_1 = (Px-Rx1)**lx1 * sp.exp(-zeta1*(Px-Rx1)**2)
    orbital_2 = (Px-Rx2)**lx2 * sp.exp(-zeta2*(Px-Rx2)**2)
    orbital_3 = sp.exp(-zeta3*(Px-Rx3)**2) # Always s-type so lx3=0.

    return sp.printing.fcode(orbital_1 * sp.diff((orbital_3 * orbital_2), Px))


def apply_seperable_substitutions(string):

    # Remove newlines, spaces, and continuation characters.
    string = re.sub("\n", "", string)
    string = re.sub(" ", "", string)
    string = re.sub("@", "", string)

    # Replace xyz positions with array references.
    string = string.replace("Px", "xyz(:)")

    # Replace angular momenta with array references.
    string = string.replace("lx1", "l1(:)")
    string = string.replace("lx2", "l2(:)")
    string = string.replace("lx3", "l3(:)")

    # Replace A,B,C site positions with array references.
    string = string.replace("Rx1", "A(:)")
    string = string.replace("Rx2", "B(:)")
    string = string.replace("Rx3", "C(:)")

    # Replace zeta with alpha (a).
    string = string.replace("zeta", "a")

    return string


def apply_paired_seperable_substitutions(string):

    # Remove newlines, spaces, and continuation characters.
    string = re.sub("\n", "", string)
    string = re.sub(" ", "", string)
    string = re.sub("@", "", string)

    # Replace xyz positions with array references.
    string = string.replace("Px", "xyz(i)")
    string = string.replace("Py", "xyz(j)")

    # Replace angular momenta with array references.
    string = string.replace("lx1", "l1(i)")
    string = string.replace("lx2", "l2(i)")
    string = string.replace("ly1", "l1(j)")
    string = string.replace("ly2", "l2(j)")

    # Replace A,B,C site positions with array references.
    string = string.replace("Rx1", "A(i)")
    string = string.replace("Rx2", "B(i)")
    string = string.replace("Rx3", "C(i)")
    string = string.replace("Ry1", "A(j)")
    string = string.replace("Ry2", "B(j)")
    string = string.replace("Ry3", "C(j)")

    # Replace zeta with alpha (a).
    string = string.replace("zeta", "a")

    return string


def apply_inseperable_substitutions(string, vectorize):

    # Remove newlines, spaces, and continuation characters.
    string = re.sub("\n", "", string)
    string = re.sub(" ", "", string)
    string = re.sub("@", "", string)

    # Replace xyz positions with array references.
    string = string.replace("Px", "xyz(1)")
    string = string.replace("Py", "xyz(2)")
    string = string.replace("Pz", "xyz(3)")

    # Replace angular momenta with array references.
    string = string.replace("lx1", "l1(1)")
    string = string.replace("ly1", "l1(2)")
    string = string.replace("lz1", "l1(3)")
    string = string.replace("lx2", "l2(1)")
    string = string.replace("ly2", "l2(2)")
    string = string.replace("lz2", "l2(3)")
    string = string.replace("lx3", "l3(1)")
    string = string.replace("ly3", "l3(2)")
    string = string.replace("lz3", "l3(3)")

    # Replace A,B,C site positions with array references.
    string = string.replace("Rx1", "A(1)")
    string = string.replace("Ry1", "A(2)")
    string = string.replace("Rz1", "A(3)")
    string = string.replace("Rx2", "B(1)")
    string = string.replace("Ry2", "B(2)")
    string = string.replace("Rz2", "B(3)")
    string = string.replace("Rx3", "C(1)")
    string = string.replace("Ry3", "C(2)")
    string = string.replace("Rz3", "C(3)")

    # Pull out the exponential terms for the second site and attach them to
    #   the first site. Take special note that this replacement uses "zeta"
    #   because we have not yet performed the zeta substitution.
    string = string.replace("*exp(-zeta2*sum((xyz(:)-B(:))**2))","")
    string = string.replace("exp(-zeta1*sum((xyz(:)-A(:))**2))",
            "exp(-zeta1*sum((xyz(:)-A(:))**2))*"+
            "exp(-zeta2*sum((xyz(:)-B(:))**2))")

    # Replace zeta with alpha (a).
    string = string.replace("zeta", "a")

    # If a request to vectorize some expressions was made, then do it.
    if (vectorize):

        # Replace sums and products with vector expressions.
        string = string.replace("((xyz(1)-A(1))**2+" +
                                "(xyz(2)-A(2))**2+" +
                                "(xyz(3)-A(3))**2)",
                                "sum((xyz(:)-A(:))**2)")
        string = string.replace("((xyz(1)-B(1))**2+" +
                                "(xyz(2)-B(2))**2+" +
                                "(xyz(3)-B(3))**2)",
                                "sum((xyz(:)-B(:))**2)")
        string = string.replace("((xyz(1)-C(1))**2+" +
                                "(xyz(2)-C(2))**2+" +
                                "(xyz(3)-C(3))**2)",
                                "sum((xyz(:)-C(:))**2)")
        string = string.replace("(xyz(1)-A(1))**l1(1)*" +
                                "(xyz(2)-A(2))**l1(2)*" +
                                "(xyz(3)-A(3))**l1(3)",
                                "product((xyz(:)-A(:))**l1(:))")
        string = string.replace("(xyz(1)-B(1))**l2(1)*" +
                                "(xyz(2)-B(2))**l2(2)*" +
                                "(xyz(3)-B(3))**l2(3)",
                                "product((xyz(:)-B(:))**l2(:))")
        string = string.replace("(xyz(1)-C(1))**l3(1)*" +
                                "(xyz(2)-C(2))**l3(2)*" +
                                "(xyz(3)-C(3))**l3(3)",
                                "product((xyz(:)-C(:))**l3(:))")

        # Condense expressions of the form:
        #   (2*xyz(#)-2*B(#))*product(xyz(:)-B(:)**l2(:))/(xyz(#)-B(#)) where
        #   "#" is a 1, 2, or 3 and the position of the term with the "2*"
        #   factor can be different in the expanded form of the product term.
        #   So, I have yet to develop a simple regular expression replacement
        #   scheme for the general case and so to move more quickly I'm just
        #   going to explicitly condense each known form. The one trick is
        #   that this condensation must be done *after* the previous
        #   condensation of the exponential terms because they get in the way
        #   of this search-and-replace. Also note, the #3 term has already
        #   been partially reduced.
        string = string.replace("(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*"+
                "(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))",
                "2*product((xyz(:)-B(:))**l2(:))")
        string = string.replace("(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*"+
                "(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))",
                "2*product((xyz(:)-B(:))**l2(:))")
        string = string.replace("product((xyz(:)-B(:))**l2(:))*"+
                "(2*xyz(3)-2*B(3))/(xyz(3)-B(3))",
                "2*product((xyz(:)-B(:))**l2(:))")

        # As above except here we condense expressions of the form:
        #   (2*xyz(#)-2*B(#))**2*product(xyz(:)-B(:)**l2(:)). Also this case,
        #   the #3 term has already been partially reduced.
        string = string.replace("(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))**2*"+
                "(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)",
                "4*product((xyz(:)-B(:))**l2(:))*(xyz(1)-B(1))**2")
        string = string.replace("(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*"+
                "(2*xyz(2)-2*B(2))**2*(xyz(3)-B(3))**l2(3)",
                "4*product((xyz(:)-B(:))**l2(:))*(xyz(2)-B(2))**2")
        string = string.replace("product((xyz(:)-B(:))**l2(:))*"+
                "(2*xyz(3)-2*B(3))**2",
                "4*product((xyz(:)-B(:))**l2(:))*(xyz(3)-B(3))**2")

        # Precompute sums and products as possible.
        string = string.replace("product((xyz(:)-A(:))**l1(:))",
                "product_xyz_A_l1")
        string = string.replace("product((xyz(:)-B(:))**l2(:))",
                "product_xyz_B_l2")
        string = string.replace("product((xyz(:)-C(:))**l3(:))",
                "product_xyz_C_l3")
        string = string.replace("sum((xyz(:)-A(:))**2)",
                "sum_xyz_A_2")
        string = string.replace("sum((xyz(:)-B(:))**2)",
                "sum_xyz_B_2")
        string = string.replace("sum((xyz(:)-C(:))**2)",
                "sum_xyz_C_2")

    return string


def main():

    # Set default values from the user resource control file.
    settings = init_environment()

    # Parse the command line.
    args = parse_command_line()

    # Reconcile the default settings and the given arguments
    settings = reconcile(args, settings)

    # Proceed with computing and printing the results.
    sp.init_printing()
    f = open("temp_Fortran_code.f90", "w")

    # Manage the two-center overlap integrals.
    if (settings.overlap):
        string = overlap()
        string = apply_seperable_substitutions(string)
        lib.print_cont_string(string, 80, 3, f, True)

    if (settings.kinetic):
        string = kinetic_energy()
        string = apply_inseperable_substitutions(string, False)
        lib.print_cont_string(string, 80, 3, f, True)

    if (settings.nuclear):
        string = nuclear()
        string = apply_inseperable_substitutions(string, False)
        lib.print_cont_string(string, 80, 3, f, True)

    if (settings.electron):
        string = electron()
        string = apply_seperable_substitutions(string)
        lib.print_cont_string(string, 80, 3, f, True)

    if (settings.momentum):
        string = momentum()
        string = apply_seperable_substitutions(string)
        lib.print_cont_string(string, 80, 3, f, True)

    if (settings.massvel):
        # Create integral formulas for 1D fourth derivatives and
        #   2D products of second derivatives.
        string1, string2 = massvel()
        string1 = apply_seperable_substitutions(string1)
        string2 = apply_paired_seperable_substitutions(string2)
        lib.print_cont_string(string1, 80, 3, f, True)
        lib.print_cont_string(string2, 80, 3, f, True)

        # Uncomment below for the "full" 3D mesh case.
        #string = massvel()
        #string = apply_inseperable_substitutions(string, False)
        #lib.print_cont_string(string, 80, 3, f, True)

    if (settings.dipole):
        string = dipole()
        string = apply_inseperable_substitutions(string, False)
        lib.print_cont_string(string, 80, 3, f, True)

    if (settings.dkinetic):
        string1, string2 = derivative_kinetic_energy()
        string1 = apply_seperable_substitutions(string1)
        string2 = apply_paired_seperable_substitutions(string2)
        lib.print_cont_string(string1, 80, 3, f, True)
        lib.print_cont_string(string2, 80, 3, f, True)

    if (settings.dnuclear):
        string1, string2, string3 = derivative_nuclear()
        string1 = apply_inseperable_substitutions(string1, False)
        string2 = apply_inseperable_substitutions(string2, False)
        string3 = apply_inseperable_substitutions(string3, False)
        lib.print_cont_string(string1, 80, 3, f, True)
        lib.print_cont_string(string2, 80, 3, f, True)
        lib.print_cont_string(string3, 80, 3, f, True)

    if (settings.delectron):
        string = derivative_electron()
        string = apply_seperable_substitutions(string)
        lib.print_cont_string(string, 80, 3, f, True)


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or a request
    #   to import information from external modules. Only now do we actually
    #   start running the program. In this case, there is no program. This
    #   should only be imported as a module.
    main()
