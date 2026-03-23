#!/usr/bin/env python3

######################################################################
#
# PROGRAM: uolcao
#
# This is the driver script for the unified OLCAO (Orthogonalized
#   Linear Combination of Atomic Orbitals) electronic structure code.
#   It orchestrates the execution of Fortran programs that perform
#   self-consistent field (SCF) calculations and a variety of
#   post-SCF property calculations. The script manages:
#     - Input file preparation and staging
#     - Selection of the correct Fortran executable (gamma vs.
#       non-gamma k-point versions)
#     - Intermediate/temporary directory management
#     - Output file collection and renaming
#     - Checkpointing (completed calculations are skipped on
#       restart)
#     - Lock file management to prevent concurrent runs
#
# The script supports a rich set of calculation types including:
#   density of states (DOS), bond order, dipole moment,
#   polarization (MTOP), optical properties (linear and
#   nonlinear), photo-absorption cross section (PACS), sigma(E),
#   symmetric band structure (SYBD), forces, field (charge
#   density, potential, wave function), and local environment
#   (LoEn) analysis.
#
# Each calculation type can be run as:
#   - An SCF-only job (e.g., -scfdos): performs SCF with the
#     property-specific basis, then computes the property
#   - A post-SCF job (e.g., -dos): performs SCF with one basis
#     and the property calculation with a (possibly different)
#     basis
#
# Three basis set sizes are available:
#   MB = Minimal Basis, FB = Full Basis, EB = Extended Basis
#
######################################################################

import argparse as ap
import os
import shutil
import subprocess
import sys
from datetime import datetime


# ------------------------------------------------------------------ #
#                     File Name Components                            #
# ------------------------------------------------------------------ #

class FileNames:
    """Container for all file name components used throughout the
    program. Centralizing them here provides a quick overview of
    every file the program touches and makes renaming easy.

    File names are built by concatenating these components. For
    example, the SCF potential file is:
        edge_ + initPot + basis + dat
        => "gs_scfV-fb.dat"

    All edge-dependent output files are prefixed with the edge
    name when actually used."""

    def __init__(self):

        # Miscellaneous control and tracking files.
        self.olcao_kill = "OLCAOkill"
        self.olcao_lock = "OLCAOlock"
        self.runtime = "runtime"
        self.intermediate = "intermediate"
        self.energy = "enrg"
        self.iteration = "iter"
        self.kp_scf = "kp-scf"
        self.kp_pscf = "kp-pscf"
        # Magnetic moment. For spin polarized calculations only.
        self.moment = "mgmo"
        self.init_pot = "scfV"
        self.structure = "structure"
        self.olcao = "olcao"
        self.atom_pos = "atomPos"
        self.lattice = "lattice"

        # Replacements (SCF vs post-SCF tags).
        self.scf = "scf"
        self.pscf = "pscf"

        # File nature components (property names).
        self.dos = "dos"
        self.bond = "bond"
        self.dimo = "dimo"
        self.optc = "optc"
        self.pacs = "pacs"
        self.nlop = "nlop"
        self.sige = "sige"
        self.sybd = "sybd"
        self.force = "force"
        self.field = "field"
        self.mtop = "mtop"
        self.loen = "loen"
        self.vdim = "vdim"

        # Auxiliary suffix components (property sub-types).
        self.loci = ".loci"
        self.tot = ".t"        # total
        self.par = ".p"        # partial
        self.bo3c = ".3c"
        self.elf = ".elf"
        self.eps1 = ".eps1"
        self.eps1i = ".eps1i"
        self.eps2 = ".eps2"
        self.chi2 = ".chi2"
        self.chi1 = ".chi1"
        self.cond = ".cond"    # Optical conductivity
        self.kkc = ".kkc"      # Partial kramers-kronig
        self.extnct = ".kext"  # Extinction coefficient
        self.refrac = ".nref"  # Refractive index
        self.absorp = ".aabs"  # Absorption coefficient
        self.reflec = ".Rref"  # Reflectivity
        self.pot = ".pot"
        self.rho = ".rho"
        self.wave = ".wave"
        self.up = ".up"
        self.dn = ".dn"
        self.up_p_dn = ".up+dn"
        self.up_m_dn = ".up-dn"
        # A suffix of a suffix. Minus the neutral atom.
        self.minus_neut = "-N"
        self.profile = ".prof"
        self.bs_comp = ".bc"

        # File name extensions.
        self.plot = ".plot"
        self.raw = ".raw"
        self.out = ".out"
        self.dat = ".dat"
        self.hdf = ".hdf5"
        self.xdmf3 = ".xdmf3"


# ------------------------------------------------------------------ #
#                   Job Type Definitions                              #
# ------------------------------------------------------------------ #

# Information about which job to run for a given jobID.
# SCF only = 0
# dosSCF=101;   dosPSCF=201
# bondSCF=102;  bondPSCF=202
# dimoSCF=103;  dimoPSCF=203
# optcSCF=104;  optcPSCF=204
# pacsSCF=105;  pacsPSCF=205
# nlopSCF=106;  nlopPSCF=206
# sigeSCF=107;  sigePSCF=207
# sybdSCF=108;  sybdPSCF=208
# forceSCF=109; forcePSCF=209
# fieldSCF=110; fieldPSCF=210
# mtopSCF=111;  mtopPSCF=211
# loenPSCF=311

# Each entry: (jobID, jobName, default_scf_basis, default_pscf_basis)
# A value of None means "do not override if already set".
JOB_DEFS = {
    # SCF-only jobs (100-series): only scf basis is defaulted.
    "scfdos":   (101, "scf+dos",   "fb", None),
    "scfbond":  (102, "scf+bond",  "mb", None),
    "scfdimo":  (103, "scf+dimo",  "fb", None),
    "scfoptc":  (104, "scf+optc",  "eb", None),
    "scfpacs":  (105, "scf+pacs",  "eb", None),
    "scfnlop":  (106, "scf+nlop",  "eb", None),
    "scfsige":  (107, "scf+sige",  "fb", None),
    "scfsybd":  (108, "scf+sybd",  "fb", None),
    "scfforce": (109, "scf+force", "fb", None),
    "scffield": (110, "scf+field", "fb", None),
    "scfmtop":  (111, "scf+mtop",  "fb", None),
    # Post-SCF jobs (200-series): both bases defaulted.
    "dos":   (201, "pscf+dos",   "fb", "fb"),
    "bond":  (202, "pscf+bond",  "fb", "mb"),
    "dimo":  (203, "pscf+dimo",  "mb", "mb"),
    "optc":  (204, "pscf+optc",  "fb", "eb"),
    "pacs":  (205, "pscf+pacs",  "fb", "eb"),
    "nlop":  (206, "pscf+nlop",  "fb", "eb"),
    "sige":  (207, "pscf+sige",  "fb", "fb"),
    "sybd":  (208, "pscf+sybd",  "fb", "fb"),
    "force": (209, "pscf+force", "fb", "fb"),
    "field": (210, "pscf+field", "fb", "fb"),
    "mtop":  (211, "pscf+mtop",  "fb", "fb"),
    # Special jobs (300-series).
    "loen":  (311, "loen",        None, None),
}


# ------------------------------------------------------------------ #
#                     Edge Code Mapping                               #
# ------------------------------------------------------------------ #

# Map edge string to (QN_n, QN_l) quantum numbers.
EDGE_MAP = {
    "gs": (0, 0),
    "1s": (1, 0),
    "2s": (2, 0), "2p": (2, 1),
    "3s": (3, 0), "3p": (3, 1), "3d": (3, 2),
    "4s": (4, 0), "4p": (4, 1), "4d": (4, 2), "4f": (4, 3),
    "5s": (5, 0), "5p": (5, 1), "5d": (5, 2),
    "6s": (6, 0), "6p": (6, 1),
    "7s": (7, 0),
}

# Map basis string to numeric code.
BASIS_CODE_MAP = {"mb": 1, "fb": 2, "eb": 3}


# ------------------------------------------------------------------ #
#                       Script Settings                               #
# ------------------------------------------------------------------ #

class ScriptSettings():
    """The instance variables of this object are the user settings
    that control the program. The variable values are pulled from
    a list that is created within a resource control file and that
    are then reconciled with command line parameters."""


    def __init__(self):
        """Define default values for the parameters by pulling
        them from the resource control file in the default
        location: $OLCAO_RC/uolcaorc.py or from the current
        working directory if a local copy of uolcaorc.py is
        present."""

        # Read default variables from the resource control file.
        rc_dir = os.getenv('OLCAO_RC')
        if not rc_dir:
            sys.exit(
                "Error: $OLCAO_RC is not set. See instructions."
            )
        sys.path.insert(1, rc_dir)
        from uolcaorc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the rc file.
        self.reconcile(args)

        # Record the command line parameters that were used.
        self.recordCLP()


    def assign_rc_defaults(self, default_rc):

        # Optical properties x,y,z serial vs. simultaneous.
        self.serialxyz = default_rc["serialxyz"]

        # Valgrind execution mode.
        self.valgrind = default_rc["valgrind"]


    def parse_command_line(self):

        # Create the parser tool.
        prog_name = "uolcao"

        description_text = """
Driver script for the unified OLCAO electronic structure code.

Orchestrates SCF and post-SCF calculations by managing input
files, selecting the correct Fortran executable, and collecting
output. Supports checkpointing: completed calculations are
skipped on restart.

EDGE is either "gs", "1s", "2s", "2p", "3s", "3p", "3d",
  "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p",
  or "7s".
BASIS is either "MB" (minimal), "FB" (full), or "EB" (extended).
"""

        epilog_text = """
PROPERTY DESCRIPTIONS:

-dos / -scfdos:
  Density of states. Default basis: FB. Outputs include TDOS
  localization index plot files (directly plottable) and a PDOS
  raw file (postprocess with makePDOS).

-bond / -scfbond:
  Bond order and Q*. Default basis: -bond uses SCF FB, post-SCF
  MB; -scfbond uses MB. Output includes bond and Q* raw files
  (postprocess with makeBOND).

-dimo / -scfdimo:
  Dipole moment. Default basis: MB. Output is a total dipole
  moment output file.

-mtop / -scfmtop:
  Polarization (multipole). Default basis: FB. Output is a
  polarization output file.

-optc / -scfoptc:
  Optical properties. Default basis: -optc uses SCF FB, post-SCF
  EB; -scfoptc uses EB. Output includes x,y,z decompositions of
  optical conductivity, epsilon1, epsilon2, energy loss function
  (ELF), refractive index, absorption coefficient, reflectivity,
  and imaginary epsilon1. Also a total eps1, eps2, and ELF file.

-pacs / -scfpacs:
  Photo-absorption cross section. EDGE is REQUIRED. Default
  basis: -pacs uses SCF FB, post-SCF EB; -scfpacs uses EB.

-nlop / -scfnlop:
  Nonlinear optical properties. Default basis: -nlop uses SCF FB,
  post-SCF EB; -scfnlop uses EB.

-sige / -scfsige:
  Sigma(E) curve of optical transitions near the Fermi energy.
  Default basis: FB.

-sybd / -scfsybd:
  Symmetric band structure. Default basis: FB. Should be done as
  a separate post-SCF calculation because it uses a k-point set
  that usually should not be used for the SCF part.

-force / -scfforce:
  Forces between atoms. Still in development; experimental.

-field / -scffield:
  Complex wave function (or real for Gamma k-point), charge
  density rho, electronic potential, spin-dependent wave function
  charge density and potential, and the difference between
  interacting and non-interacting system wave function, rho, and
  potential (with spin dependency). Also computes 1-D axial
  profiles and charge centers. Default basis: FB.

-loen:
  Local environment analysis. Quantifies the local environment
  around each atom via a cutoff factor and bispectrum components.

CHECKPOINTING:
  Checkpointing is used internally for SCF and post-SCF
  calculations. Any completed checkpointed calculations will be
  skipped if found. For example, if a job dies part way through
  SCF integrals, the job will not recompute them on restart.
  Similarly, if an SCF calculation for a specific basis is already
  complete, requesting another calculation that requires the same
  basis SCF will skip it.

Please contact Paul Rulis (rulisp@umkc.edu) regarding questions.
Defaults are given in ./uolcaorc.py or $OLCAO_RC/uolcaorc.py.
"""

        parser = ap.ArgumentParser(
            prog=prog_name,
            formatter_class=ap.RawDescriptionHelpFormatter,
            description=description_text,
            epilog=epilog_text,
        )

        # Add arguments to the parser.
        self.add_parser_arguments(parser)

        # Parse the arguments and return the results.
        return parser.parse_args()


    def add_parser_arguments(self, parser):

        # SCF basis selection.
        parser.add_argument(
            '-scf', dest='scf_basis', type=str,
            default=None, metavar='BASIS',
            help='Basis for the SCF portion: MB, FB, EB, or '
                 'NO (skip SCF). Default depends on job type.'
        )

        # Post-SCF basis selection.
        parser.add_argument(
            '-pscf', dest='pscf_basis', type=str,
            default=None, metavar='BASIS',
            help='Basis for the post-SCF portion: MB, FB, EB. '
                 'Default depends on job type.'
        )

        # Serial xyz treatment for optical properties.
        parser.add_argument(
            '-serialxyz', dest='serialxyz',
            action='store_true', default=self.serialxyz,
            help='Compute optical x,y,z components in serial '
                 'fashion (conserves memory, costs time). '
                 f'Default: {self.serialxyz}'
        )

        # Valgrind execution mode.
        parser.add_argument(
            '-valgrind', dest='valgrind',
            action='store_true', default=self.valgrind,
            help='Run Fortran programs under valgrind with '
                 '--leak-check. '
                 f'Default: {self.valgrind}'
        )

        # Create a mutually exclusive group for job types.
        #   Each flag stores two pieces of information:
        #     1. Which calculation was requested (identified
        #        by which dest is non-None after parsing).
        #     2. Which edge to use (the value stored at that
        #        dest; defaults to "gs" if not supplied).
        #   For most flags, the edge is optional (nargs='?',
        #   const='gs'). For -pacs and -scfpacs the edge is
        #   required (no nargs, so argparse demands a value).
        jobs = parser.add_mutually_exclusive_group()

        # --- SCF-only job types (100-series) --- #

        jobs.add_argument(
            '-scfdos', dest='scfdos', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='SCF + DOS calculation. '
                 'Default SCF basis: FB.'
        )
        jobs.add_argument(
            '-scfbond', dest='scfbond', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='SCF + bond order/Q*. '
                 'Default SCF basis: MB.'
        )
        jobs.add_argument(
            '-scfdimo', dest='scfdimo', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='SCF + dipole moment. '
                 'Default SCF basis: FB.'
        )
        jobs.add_argument(
            '-scfoptc', dest='scfoptc', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='SCF + optical properties. '
                 'Default SCF basis: EB.'
        )
        jobs.add_argument(
            '-scfpacs', dest='scfpacs', type=str,
            default=None, metavar='EDGE',
            help='SCF + PACS. EDGE is required. '
                 'Default SCF basis: EB.'
        )
        jobs.add_argument(
            '-scfnlop', dest='scfnlop', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='SCF + nonlinear optical properties. '
                 'Default SCF basis: EB.'
        )
        jobs.add_argument(
            '-scfsige', dest='scfsige', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='SCF + sigma(E). '
                 'Default SCF basis: FB.'
        )
        jobs.add_argument(
            '-scfsybd', dest='scfsybd', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='SCF + symmetric band structure. '
                 'Default SCF basis: FB.'
        )
        jobs.add_argument(
            '-scfforce', dest='scfforce', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='SCF + forces (experimental). '
                 'Default SCF basis: FB.'
        )
        jobs.add_argument(
            '-scffield', dest='scffield', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='SCF + field. '
                 'Default SCF basis: FB.'
        )
        jobs.add_argument(
            '-scfmtop', dest='scfmtop', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='SCF + polarization. '
                 'Default SCF basis: FB.'
        )

        # --- Post-SCF job types (200-series) --- #

        jobs.add_argument(
            '-dos', dest='dos', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='Post-SCF DOS. '
                 'Default SCF basis: FB, post-SCF: FB.'
        )
        jobs.add_argument(
            '-bond', dest='bond', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='Post-SCF bond order/Q*. '
                 'Default SCF basis: FB, post-SCF: MB.'
        )
        jobs.add_argument(
            '-dimo', dest='dimo', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='Post-SCF dipole moment. '
                 'Default SCF: MB, post-SCF: MB.'
        )
        jobs.add_argument(
            '-optc', dest='optc', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='Post-SCF optical properties. '
                 'Default SCF: FB, post-SCF: EB.'
        )
        jobs.add_argument(
            '-pacs', dest='pacs', type=str,
            default=None, metavar='EDGE',
            help='Post-SCF PACS. EDGE is required. '
                 'Default SCF: FB, post-SCF: EB.'
        )
        jobs.add_argument(
            '-nlop', dest='nlop', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='Post-SCF nonlinear optical. '
                 'Default SCF: FB, post-SCF: EB.'
        )
        jobs.add_argument(
            '-sige', dest='sige', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='Post-SCF sigma(E). '
                 'Default SCF: FB, post-SCF: FB.'
        )
        jobs.add_argument(
            '-sybd', dest='sybd', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='Post-SCF symmetric band structure. '
                 'Default SCF: FB, post-SCF: FB.'
        )
        jobs.add_argument(
            '-force', dest='force', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='Post-SCF forces (experimental). '
                 'Default SCF: FB, post-SCF: FB.'
        )
        jobs.add_argument(
            '-field', dest='field', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='Post-SCF field. '
                 'Default SCF: FB, post-SCF: FB.'
        )
        jobs.add_argument(
            '-mtop', dest='mtop', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='Post-SCF polarization. '
                 'Default SCF: FB, post-SCF: FB.'
        )
        jobs.add_argument(
            '-loen', dest='loen', nargs='?',
            const='gs', default=None, metavar='EDGE',
            help='Local environment analysis.'
        )


    def reconcile(self, args):
        """Reconcile command line arguments with rc defaults.

        This method resolves the job type, edge, and basis set
        selections. The logic mirrors the original Perl parseCLP
        subroutine where each job type flag sets two variables:
          1. The calculation type (jobID / jobName), identified
             by which argument dest is non-None.
          2. The edge (the value stored at that dest).
        Default bases are applied only if the user did not set
        them explicitly via -scf or -pscf."""

        self.serialxyz = args.serialxyz
        self.valgrind = args.valgrind

        # Track whether the user explicitly set each basis.
        scf_explicit = args.scf_basis is not None
        pscf_explicit = args.pscf_basis is not None

        # Start with the user-supplied basis or defaults.
        # String for which scf basis to use (default full).
        self.basis_scf = (
            args.scf_basis.lower()
            if scf_explicit else "fb"
        )
        # String for which pscf basis to use (default none).
        self.basis_pscf = (
            args.pscf_basis.lower()
            if pscf_explicit else "no"
        )

        # Scan the job-type dests to find which one was set.
        #   Each job-type flag has its own dest (e.g.,
        #   args.dos, args.scfdos). Exactly one will be
        #   non-None if a job flag was given, or all will be
        #   None for a default SCF-only run.
        job_type = None
        edge = None
        for name in JOB_DEFS:
            val = getattr(args, name, None)
            if val is not None:
                job_type = name
                edge = val
                break

        if job_type is None:
            # No job type specified: default SCF-only.
            self.job_id = 0
            self.job_name = "scf"
            self.edge = "gs"
        else:
            job_id, job_name, def_scf, def_pscf = (
                JOB_DEFS[job_type]
            )
            self.job_id = job_id
            self.job_name = job_name
            self.edge = edge

            # Apply default bases only if the user did not
            #   set them explicitly on the command line.
            if not scf_explicit and def_scf is not None:
                self.basis_scf = def_scf
            if not pscf_explicit and def_pscf is not None:
                self.basis_pscf = def_pscf

        # Determine numeric basis codes.
        self.basis_code_scf = BASIS_CODE_MAP.get(
            self.basis_scf, 0
        )
        self.basis_code_pscf = BASIS_CODE_MAP.get(
            self.basis_pscf, 0
        )

        # Determine edge quantum numbers.
        if self.edge not in EDGE_MAP:
            sys.exit(
                f"Error: unknown edge '{self.edge}'. "
                f"Valid: "
                f"{', '.join(EDGE_MAP.keys())}"
            )
        self.qn_n, self.qn_l = EDGE_MAP[self.edge]


    def recordCLP(self):
        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime(
                "%b. %d, %Y: %H:%M:%S"
            )
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write("Cmnd:")
            for argument in sys.argv:
                cmd.write(f" {argument}")
            cmd.write("\n\n")


# ------------------------------------------------------------------ #
#                      File Operation Helpers                         #
# ------------------------------------------------------------------ #

def safe_copy(src, dst, runtime_fh):
    """Copy src to dst. If src does not exist but dst does,
    log and skip. Dies on failure."""
    if not os.path.exists(src) and os.path.exists(dst):
        runtime_fh.write(
            f"{src} not found but {dst} exists. "
            f"No copy done.\n"
        )
        return
    try:
        shutil.copy2(src, dst)
    except Exception as e:
        sys.exit(f"Cannot copy {src} to {dst}. {e}")


def safe_move(src, dst, runtime_fh):
    """Move src to dst. If src does not exist but dst does,
    log and skip. Dies on failure."""
    if not os.path.exists(src) and os.path.exists(dst):
        runtime_fh.write(
            f"{src} not found but {dst} exists. "
            f"No move done.\n"
        )
        return
    try:
        shutil.move(src, dst)
    except Exception as e:
        sys.exit(f"Cannot move {src} to {dst}. {e}")


def safe_append(src, dst, skip_lines, runtime_fh):
    """Append src to dst, skipping the first skip_lines lines
    of src. If dst does not exist, copy src entirely. If src
    does not exist but dst does, log and skip."""
    if not os.path.exists(dst):
        try:
            shutil.copy2(src, dst)
        except Exception as e:
            sys.exit(
                f"Cannot full append {src} to {dst}. {e}"
            )
        return
    if not os.path.exists(src) and os.path.exists(dst):
        runtime_fh.write(
            f"{src} not found but {dst} exists. "
            f"No append done.\n"
        )
        return
    try:
        with open(src, 'r') as f_in:
            lines = f_in.readlines()
        # skip_lines is 1-based (like tail -n +N).
        with open(dst, 'a') as f_out:
            for line in lines[skip_lines - 1:]:
                f_out.write(line)
    except Exception as e:
        sys.exit(
            f"Cannot append {src} to {dst}. {e}"
        )


def safe_delete(path):
    """Delete a file. Dies on failure."""
    try:
        os.remove(path)
    except Exception as e:
        sys.exit(f"Cannot delete {path}. {e}")


def olcao_exit(msg, runtime_fh):
    """Record a reason for exiting to the RUNTIME file and
    exit the script."""
    runtime_fh.write(f"{msg}\n")
    runtime_fh.close()
    sys.exit(msg)


# ------------------------------------------------------------------ #
#                   XC Code / Spin Polarization                       #
# ------------------------------------------------------------------ #

def xc_code_spin(proj_home, inputs, fn):
    """Determine the XC code from olcao.dat and look up the
    spin polarization value from the XC functional database.

    This subroutine reads the XC_CODE from the olcao.dat input
    file and cross-references it with the xc_code.dat database
    to determine whether the calculation is spin polarized.

    Returns:
        int: 1 for non-spin-polarized, 2 for spin-polarized.
    """

    # Locate and open olcao.dat for reading.
    olcao_path = os.path.join(
        proj_home, f"{fn.olcao}{fn.dat}"
    )
    if not os.path.exists(olcao_path):
        olcao_path = os.path.join(
            inputs, f"{fn.olcao}{fn.dat}"
        )
    if not os.path.exists(olcao_path):
        sys.exit(
            f"No {fn.olcao}{fn.dat} input file found "
            f"in {proj_home} or {inputs}."
        )

    # Read the OLCAODAT file looking for XC_CODE.
    xc_code = None
    with open(olcao_path, 'r') as f:
        for line in f:
            if "XC_CODE" in line:
                next_line = next(f).strip()
                xc_code = next_line.split()[0]
                break

    if xc_code is None:
        sys.exit("XC_CODE not found in olcao.dat.")

    # Open XC functional database file for reading.
    olcao_data = os.environ.get('OLCAO_DATA', '')
    db_path = os.path.join(olcao_data, "xc_code.dat")
    if not os.path.exists(db_path):
        sys.exit(
            "Could not open XC functional database file "
            f"for reading: {db_path}"
        )

    # Search the database for the matching XC code.
    spin_value = 1  # Default non-spin-polarized.
    with open(db_path, 'r') as f:
        for line in f:
            if xc_code in line:
                parts = line.split()
                spin_value = int(parts[2])
                break

    return spin_value


# ------------------------------------------------------------------ #
#                  Temporary Directory Setup                           #
# ------------------------------------------------------------------ #

def get_temp_dir():
    """Determine the location for the temporary storage
    directory based on the name of the current directory.

    The function finds the user's home directory component in
    the current path and appends everything after it to the
    $OLCAO_TEMP base directory. This creates a parallel
    directory structure under the temp location."""
    directory = os.getcwd()
    directories = directory.split('/')
    user = os.environ.get('USER', '')
    olcao_temp = os.environ.get('OLCAO_TEMP', '')

    # Look for the directory where the username is specified.
    #   Once found, append the remaining directories to the
    #   temp directory.
    found = False
    temp = olcao_temp
    for d in directories:
        if found:
            temp = os.path.join(temp, d)
        if d == user:
            found = True

    return temp


def init_directories(fn):
    """Initialize the directories used by the script.

    Sets up the project home directory, inputs directory, and
    temporary storage directory. Creates a symbolic link named
    'intermediate' pointing to the temp directory. If an
    existing link points to the wrong location, it is renamed
    with a FIXME suffix and a new link is created.

    Returns:
        tuple: (proj_home, inputs, temp) directory paths.
    """

    # Define the project's home directory.
    proj_home = os.getcwd()

    # Define the location where the input files are kept.
    inputs = os.path.join(proj_home, "inputs")

    # Determine the location for the temporary storage
    #   directory based on the name of the current directory.
    temp = get_temp_dir()

    # Create the temp directory if it does not already exist.
    os.makedirs(temp, exist_ok=True)

    # If a link to the temp directory does not yet exist, then
    #   create one. If a link already exists but points to the
    #   wrong place, then rename the existing link and create
    #   a new one.
    intermediate = fn.intermediate
    if not os.path.exists(intermediate):
        os.symlink(temp, intermediate)
    else:
        # Check where the existing link points.
        if os.path.islink(intermediate):
            current_target = os.readlink(intermediate)
        else:
            current_target = None
        if current_target != temp:
            safe_move(
                intermediate,
                intermediate + "FIXME",
                sys.stderr,
            )
            os.symlink(temp, intermediate)

    return proj_home, inputs, temp


# ------------------------------------------------------------------ #
#                        IO Initialization                            #
# ------------------------------------------------------------------ #

def init_io(proj_home, fn):
    """Initialize IO for storing results to the RUNTIME file.

    Opens (or creates) the runtime log file and redirects
    stderr to it.

    Returns:
        file: The open RUNTIME file handle for writing.
    """
    runtime_path = os.path.join(proj_home, fn.runtime)

    # Open the RUNTIME file for appending (or create it).
    runtime_fh = open(runtime_path, 'a')

    return runtime_fh


# ------------------------------------------------------------------ #
#                   Executable Initialization                         #
# ------------------------------------------------------------------ #

def check_gamma_kp(kp_file, inputs):
    """Open a kpoint file and check if it requests a gamma
    kpoint (single k-point at 0,0,0).

    The kpoint file format is:
        Line 1: kpoint style code tag
        Line 2: kpoint style code value
        Line 3: integration code tag
        Line 4: integration code value
        Line 5: number of kpoints tag
        Line 6: number of kpoints value
        Line 7: explicit kpoint list header tag
        Line 8+: kpoint coordinates

    Returns:
        bool: True if this is a gamma-point-only calculation.
    """
    kp_path = os.path.join(inputs, kp_file)
    if not os.path.exists(kp_path):
        sys.exit(
            f"Cannot open {kp_path} for reading."
        )

    with open(kp_path, 'r') as f:
        lines = f.readlines()

    # Read past the kpoint style code tag and value, the
    #   integration code tag and value, and the number of
    #   kpoints tag (5 lines). Then read the number of
    #   kpoints (1 for Gamma).
    num_kp_line = lines[5].split()
    num_kp = int(num_kp_line[0])

    # Read past the explicit kpoint list header tag and
    #   then read the specific coordinates for the first
    #   kpoint.
    kp_line = lines[7].split()

    # If the number of kpoints equals 1 and the only kpoint
    #   is at 0,0,0, then we are doing a Gamma kpoint
    #   calculation.
    if (num_kp == 1
            and float(kp_line[2]) == 0
            and float(kp_line[3]) == 0
            and float(kp_line[4]) == 0):
        return True
    else:
        return False


def init_exes(settings, fn, inputs):
    """Initialize the executables and execution method.

    Determines whether to use the gamma-point version
    (guOLCAO) or the general k-point version (uOLCAO) of
    the executable based on the kpoint files.

    Returns:
        tuple: (executable, exe_mechanism, sub_mechanism)
    """
    bin_dir = os.environ.get('OLCAO_BIN', '')

    # Open each kpoint file and check if it requests a
    #   gamma kpoint.
    do_gamma_scf = check_gamma_kp(
        f"{fn.kp_scf}{fn.dat}", inputs
    )
    do_gamma_pscf = check_gamma_kp(
        f"{fn.kp_pscf}{fn.dat}", inputs
    )

    # Determine if the calculation is going to use the gamma
    #   kpoint and assign the executable accordingly.
    if do_gamma_scf and do_gamma_pscf:
        executable = "guOLCAO"
    elif not do_gamma_scf and not do_gamma_pscf:
        executable = "uOLCAO"
    else:
        # NEED TO FIX! If someone runs an SCF gamma and then
        #   a postSCF non-gamma, we need to switch executables
        #   when going from SCF to postSCF. I.e., two
        #   different executions would need to be run. For now,
        #   we just tell the user to do two separate
        #   calculations explicitly.
        sys.exit(
            "Attempting mixed gamma/non-gamma calculations. "
            "Run separately."
        )

    # Set the command line execution mechanism.
    if not settings.valgrind:
        exe_mechanism = "time"
    else:
        exe_mechanism = "time valgrind --leak-check=yes"

    # Determine the mechanism for submitting a serial OLCAO
    #   job.
    sub_mechanism = os.environ.get('OLCAO_EXE', '')

    return executable, exe_mechanism, sub_mechanism


# ------------------------------------------------------------------ #
#                     Input File Management                           #
# ------------------------------------------------------------------ #

def manage_input(settings, fn, proj_home, inputs, temp,
                 runtime_fh):
    """Check for the input files needed for this job.

    Checks for the primary 'olcao' input file and the system
    structure in the home directory; if not found, copies them
    from the inputs directory. Copies the universally required
    input files to fort.* unit numbers in the intermediate
    directory. Handles kpoint files and potential coefficient
    files based on the job type."""

    dat = fn.dat

    # Check for the primary "olcao" input file and the system
    #   structure in the home directory and if they are not
    #   there, get them from the inputs directory.
    olcao_dat = f"{fn.olcao}{dat}"
    struct_dat = f"{fn.structure}{dat}"
    if not os.path.exists(
            os.path.join(proj_home, olcao_dat)):
        safe_copy(
            os.path.join(inputs, olcao_dat),
            os.path.join(proj_home, olcao_dat),
            runtime_fh,
        )
    if not os.path.exists(
            os.path.join(proj_home, struct_dat)):
        safe_copy(
            os.path.join(inputs, struct_dat),
            os.path.join(proj_home, struct_dat),
            runtime_fh,
        )

    # Copy the universally required input files into the
    #   intermediate dir.
    safe_copy(
        os.path.join(proj_home, olcao_dat),
        "fort.5", runtime_fh,
    )
    safe_copy(
        os.path.join(proj_home, struct_dat),
        "fort.4", runtime_fh,
    )

    # Check for the set of kpoints needed for the job.
    kp_scf_dat = f"{fn.kp_scf}{dat}"
    kp_pscf_dat = f"{fn.kp_pscf}{dat}"

    if 0 <= settings.job_id < 200:
        # Doing an SCF job only.
        if not os.path.exists(
                os.path.join(proj_home, kp_scf_dat)):
            safe_copy(
                os.path.join(inputs, kp_scf_dat),
                os.path.join(proj_home, kp_scf_dat),
                runtime_fh,
            )
        safe_copy(
            os.path.join(proj_home, kp_scf_dat),
            "fort.15", runtime_fh,
        )
        get_pot_coeffs(
            settings.basis_scf, settings.edge, fn,
            proj_home, inputs, runtime_fh,
        )
    elif 200 <= settings.job_id < 300:
        # Doing a postSCF job.
        if not os.path.exists(
                os.path.join(proj_home, kp_pscf_dat)):
            safe_copy(
                os.path.join(inputs, kp_pscf_dat),
                os.path.join(proj_home, kp_pscf_dat),
                runtime_fh,
            )
        safe_copy(
            os.path.join(proj_home, kp_pscf_dat),
            "fort.16", runtime_fh,
        )
        get_pot_coeffs(
            settings.basis_pscf, settings.edge, fn,
            proj_home, inputs, runtime_fh,
        )

        # If the SCF part does not have an explicit "NO" for
        #   the basis, then we should also copy the scf
        #   kpoints.
        if settings.basis_code_scf > 0:
            if not os.path.exists(
                    os.path.join(proj_home, kp_scf_dat)):
                safe_copy(
                    os.path.join(inputs, kp_scf_dat),
                    os.path.join(proj_home, kp_scf_dat),
                    runtime_fh,
                )
            safe_copy(
                os.path.join(proj_home, kp_scf_dat),
                "fort.15", runtime_fh,
            )
    elif settings.job_id == 311:
        # A local environment (LoEn) job. Like a post SCF job.
        if not os.path.exists(
                os.path.join(proj_home, kp_pscf_dat)):
            safe_copy(
                os.path.join(inputs, kp_pscf_dat),
                os.path.join(proj_home, kp_pscf_dat),
                runtime_fh,
            )
        safe_copy(
            os.path.join(proj_home, kp_pscf_dat),
            "fort.16", runtime_fh,
        )
        get_pot_coeffs(
            settings.basis_pscf, settings.edge, fn,
            proj_home, inputs, runtime_fh,
        )

        # If the SCF part does not have an explicit "NO" for
        #   the basis, then we should also copy the scf
        #   kpoints.
        if settings.basis_code_scf > 0:
            if not os.path.exists(
                    os.path.join(proj_home, kp_scf_dat)):
                safe_copy(
                    os.path.join(inputs, kp_scf_dat),
                    os.path.join(proj_home, kp_scf_dat),
                    runtime_fh,
                )
            safe_copy(
                os.path.join(proj_home, kp_scf_dat),
                "fort.15", runtime_fh,
            )

    # Check for a potential coefficient file for SCF
    #   calculations.
    if settings.basis_code_scf > 0:
        init_pot_dat = f"{fn.init_pot}{dat}"
        if not os.path.exists(
                os.path.join(proj_home, init_pot_dat)):
            safe_copy(
                os.path.join(inputs, init_pot_dat),
                os.path.join(proj_home, init_pot_dat),
                runtime_fh,
            )


def get_pot_coeffs(task_basis, edge, fn, proj_home,
                   inputs, runtime_fh):
    """Search for a suitable set of potential coefficients.

    The search order is:
      1. Same edge and basis as the current task
      2. Same edge, alternative bases (EB, FB, MB)
      3. Ground state with the current basis
      4. Ground state, alternative bases (EB, FB, MB)
      5. Initial non-SCF potential in project home
      6. Initial non-SCF potential in inputs directory

    The first match found is copied to fort.8 for use by
    the Fortran programs."""

    dat = fn.dat
    init_pot = fn.init_pot
    edge_ = f"{edge}_"
    basis = f"-{task_basis}"
    alt_basis_list = ["-eb", "-fb", "-mb"]

    # Check for the potential from the same edge and basis
    #   of the current task.
    pot_path = os.path.join(
        proj_home, f"{edge_}{init_pot}{basis}{dat}"
    )
    if os.path.exists(pot_path):
        safe_copy(pot_path, "fort.8", runtime_fh)
        runtime_fh.write(
            f"Using {edge_}{init_pot}{basis}{dat}\n"
        )
        return

    # Check each of the alternative basis sets for the
    #   current edge.
    for alt_basis in alt_basis_list:
        pot_path = os.path.join(
            proj_home,
            f"{edge_}{init_pot}{alt_basis}{dat}",
        )
        if os.path.exists(pot_path):
            safe_copy(pot_path, "fort.8", runtime_fh)
            runtime_fh.write(
                f"Using "
                f"{edge_}{init_pot}{alt_basis}{dat}\n"
            )
            return

    # Check the ground state of the current basis.
    pot_path = os.path.join(
        proj_home, f"gs_{init_pot}{basis}{dat}"
    )
    if os.path.exists(pot_path):
        safe_copy(pot_path, "fort.8", runtime_fh)
        runtime_fh.write(
            f"Using gs_{init_pot}{basis}{dat}\n"
        )
        return

    # Check each of the alternative basis sets for the
    #   ground state.
    for alt_basis in alt_basis_list:
        pot_path = os.path.join(
            proj_home,
            f"gs_{init_pot}{alt_basis}{dat}",
        )
        if os.path.exists(pot_path):
            safe_copy(pot_path, "fort.8", runtime_fh)
            runtime_fh.write(
                f"Using gs_{init_pot}{alt_basis}{dat}\n"
            )
            return

    # Check for the initial non-scf potential in the project
    #   home directory.
    pot_path = os.path.join(
        proj_home, f"{init_pot}{dat}"
    )
    if os.path.exists(pot_path):
        runtime_fh.write(
            f"Using {init_pot}{dat} from {proj_home}.\n"
        )
        safe_copy(pot_path, "fort.8", runtime_fh)
        return

    # Check for the initial non-scf potential in the inputs
    #   directory.
    pot_path = os.path.join(inputs, f"{init_pot}{dat}")
    if os.path.exists(pot_path):
        runtime_fh.write(
            f"Using {init_pot}{dat} from {inputs}.\n"
        )
        safe_copy(pot_path, "fort.8", runtime_fh)
        return


# ------------------------------------------------------------------ #
#                     Program Execution                               #
# ------------------------------------------------------------------ #

def execute_program(job_clp, settings, fn, bin_dir,
                    executable, exe_mechanism,
                    sub_mechanism, spin_pol, runtime_fh):
    """Execute the requested Fortran program.

    Constructs and runs the command line for the OLCAO Fortran
    executable. After execution, handles secondary jobs that
    must run immediately afterwards:
      - SYBD (jobID % 100 == 8): runs makeSYBD
      - OPTC/NLOP (jobID % 100 == 4 or 6): runs OLCAOkkc
        for Kramers-Kronig conversion, and processPOPTC for
        partial optical properties if applicable.

    Finally, checks for the fort.2 success file that signals
    the Fortran executable completed without abortive error.
    """

    # If we want to run symmetric band structure calculations,
    #   then we cannot use the gamma k-point version of the
    #   program.
    exe = executable
    if (exe.startswith('g')
            and (settings.job_id == 108
                 or settings.job_id == 208)):
        exe = exe[1:]

    # Call the executable.
    cmd = (
        f"{sub_mechanism} {exe_mechanism} "
        f"{bin_dir}/{exe} {job_clp}"
    )
    result = subprocess.run(
        cmd, shell=True,
        capture_output=True, text=True,
    )
    runtime_fh.write(result.stdout)
    if result.stderr:
        runtime_fh.write(result.stderr)

    # In certain cases secondary jobs need to be run
    #   immediately afterwards.
    if settings.job_id % 100 == 8:  # SYBD
        subprocess.run(
            f"{bin_dir}/makeSYBD -dat fort.5 -out fort.20 "
            f"-raw fort.31 -plot fort.41",
            shell=True,
        )
        if spin_pol == 2:  # Spin polarized case
            subprocess.run(
                f"{bin_dir}/makeSYBD -dat fort.5 "
                f"-out fort.20 -raw fort.32 -plot fort.42",
                shell=True,
            )
    elif settings.job_id % 100 in (4, 6):  # OPTC or NLOP
        # Perform the Kramers-Kronig conversion of eps2 to
        #   eps1 for the linear optical properties and of
        #   chi2 to chi1 for the nonlinear optical properties.
        #   Also compute the energy loss function (all done in
        #   OLCAOkkc). We need to pass the program the number
        #   of lines in the input file and which set of file
        #   numbers to use. A 1 means use the spin up or
        #   default file numbers, and a 2 means use the spin
        #   down file numbers.
        result = subprocess.run(
            "wc -l fort.50",
            shell=True, capture_output=True, text=True,
        )
        optc_lines = result.stdout.split()[0]

        # The optc, poptc, and nlop calculations all need to
        #   run OLCAOkkc on the "standard" total eps2.

        # We always use the #1 file numbers.
        subprocess.run(
            f"{bin_dir}/OLCAOkkc "
            f"{optc_lines} 1 0 1",
            shell=True, capture_output=True, text=True,
        )

        # We only use the #2 set of file numbers for the
        #   spin polarized calculations.
        if spin_pol == 2:
            subprocess.run(
                f"{bin_dir}/OLCAOkkc "
                f"{optc_lines} 2 0 1",
                shell=True,
                capture_output=True, text=True,
            )

        # If we are doing the partial optical properties,
        #   call processPOPTC for the OLCAOkkc calculation.
        #   (The fort.209 file contains the information
        #   needed to produce partial epsilon1 using
        #   Kramers-Kronig conversion. It is only produced
        #   by partial optical properties calculations.)
        #   The processPOPTC will repeatedly call the
        #   OLCAOkkc program.
        if os.path.exists("fort.209"):
            if spin_pol == 1:
                subprocess.run(
                    f"{bin_dir}/processPOPTC 1",
                    shell=True,
                )
            else:
                subprocess.run(
                    f"{bin_dir}/processPOPTC 1",
                    shell=True,
                )
                subprocess.run(
                    f"{bin_dir}/processPOPTC 2",
                    shell=True,
                )

    # Check for the existence of the fort.2 file that signals
    #   completion of the fortran executable without abortive
    #   error.
    if not os.path.exists("fort.2"):
        olcao_exit(
            "Fortran success file missing. Exiting Script.",
            runtime_fh,
        )
    else:
        os.remove("fort.2")


# ------------------------------------------------------------------ #
#                     Output File Management                          #
# ------------------------------------------------------------------ #

def manage_output(settings, fn, proj_home, spin_pol,
                  runtime_fh):
    """Manage the output files produced by the Fortran programs.

    Copies and moves fort.* files to their proper named
    locations in both the intermediate directory (for
    archival) and the project home directory (for the user).
    The naming convention encodes the edge, job name, basis,
    and property type.

    The operations are organized by job type (jobID % 100)
    and handle both spin-polarized and non-spin-polarized
    cases."""

    job_id = settings.job_id
    edge_ = f"{settings.edge}_"

    # Define the basis tag.
    if 0 <= job_id < 200:
        basis = f"-{settings.basis_scf}"
    elif 200 <= job_id < 300:
        basis = f"-{settings.basis_pscf}"
    else:
        basis = "-fb"

    jn = settings.job_name
    ph = proj_home
    rt = runtime_fh

    # --- SCF output (always produced unless basis is NO) --- #
    if job_id < 200 or settings.basis_scf != "no":
        safe_append(
            "fort.7",
            f"{edge_}{jn}{basis}.{fn.iteration}.7",
            2, rt,
        )
        safe_append(
            "fort.7",
            f"{ph}/{edge_}{fn.iteration}{basis}{fn.dat}",
            2, rt,
        )
        safe_delete("fort.7")
        safe_copy(
            "fort.14",
            f"{edge_}{jn}{basis}.{fn.energy}.14",
            rt,
        )
        safe_move(
            "fort.14",
            f"{ph}/{edge_}{fn.energy}{basis}{fn.dat}",
            rt,
        )
        safe_move(
            "fort.15",
            f"{fn.kp_scf}.15",
            rt,
        )
        safe_copy(
            "fort.8",
            f"{edge_}{jn}{basis}.SCFV.8",
            rt,
        )
        safe_move(
            "fort.8",
            f"{ph}/{edge_}{fn.init_pot}{basis}{fn.dat}",
            rt,
        )
        if os.path.exists("fort.1000"):
            safe_copy(
                "fort.1000",
                f"{edge_}{jn}{basis}.iterTDOS.1000",
                rt,
            )
            safe_move(
                "fort.1000",
                f"{ph}/{edge_}{jn}{basis}"
                f".iterTDOS{fn.plot}",
                rt,
            )

    # --- Property-specific output --- #

    if job_id % 100 == 1:  # Density of states
        _manage_dos_output(
            edge_, basis, fn, ph, spin_pol, rt
        )
    elif job_id % 100 == 2:  # Bond order + Q*
        _manage_bond_output(
            edge_, basis, fn, ph, spin_pol, rt
        )
    elif job_id % 100 == 3:  # Dipole moment
        _manage_dimo_output(
            edge_, basis, fn, ph, spin_pol, rt
        )
    elif job_id % 100 == 4:  # Optical properties
        _manage_optc_output(
            edge_, basis, fn, ph, spin_pol, rt
        )
    elif job_id % 100 == 5:  # PACS
        _manage_pacs_output(
            edge_, basis, fn, ph, spin_pol, rt
        )
    elif job_id % 100 == 6:  # Nonlinear optical
        _manage_nlop_output(
            edge_, basis, fn, ph, spin_pol, rt
        )
    elif job_id % 100 == 7:  # Sigma(E)
        _manage_sige_output(
            edge_, basis, fn, ph, spin_pol, rt
        )
    elif job_id % 100 == 8:  # Symmetric band structure
        _manage_sybd_output(
            edge_, basis, fn, ph, spin_pol, rt
        )
    elif job_id % 100 == 9:  # Force
        _manage_force_output(
            edge_, basis, fn, ph, spin_pol, rt
        )
    elif job_id % 100 == 10:  # Field
        _manage_field_output(
            edge_, basis, fn, ph, spin_pol, rt
        )
    elif job_id % 100 == 11:  # Polarization (mtop)
        _manage_mtop_output(
            edge_, basis, fn, ph, spin_pol, rt
        )

    if job_id == 311:  # Local environment
        _manage_loen_output(
            edge_, basis, fn, ph, rt
        )

    # If a post SCF job was done, then do this.
    if 200 <= job_id < 300:
        safe_move("fort.16", f"{fn.kp_pscf}.16", rt)

    # If a local environment job was done, then do this.
    if job_id == 311:
        safe_move("fort.16", f"{fn.kp_pscf}.16", rt)

    # Perform these moves for all tasks.
    safe_move(
        "fort.5",
        f"{edge_}{jn}{basis}{fn.dat}.5",
        rt,
    )
    safe_move(
        "fort.4",
        f"{fn.structure}{fn.dat}.4",
        rt,
    )
    safe_copy(
        "fort.20",
        f"{edge_}{jn}{basis}{fn.out}.20",
        rt,
    )
    safe_move(
        "fort.20",
        f"{ph}/{edge_}{jn}{basis}{fn.out}",
        rt,
    )


# ---- Per-property output management helpers ---- #

def _manage_dos_output(e_, b, fn, ph, sp, rt):
    """Density of states output files."""
    if sp == 2:
        safe_copy("fort.60",
                  f"{e_}{fn.dos}{b}{fn.tot}{fn.plot}.60", rt)
        safe_move("fort.60",
                  f"{ph}/{e_}{fn.dos}{b}{fn.tot}{fn.up}"
                  f"{fn.plot}", rt)
        safe_copy("fort.70",
                  f"{e_}{fn.dos}{b}{fn.par}{fn.raw}.70", rt)
        safe_move("fort.70",
                  f"{ph}/{e_}{fn.dos}{b}{fn.par}{fn.up}"
                  f"{fn.raw}", rt)
        safe_copy("fort.80",
                  f"{e_}{fn.dos}{b}{fn.loci}{fn.plot}.80",
                  rt)
        safe_move("fort.80",
                  f"{ph}/{e_}{fn.dos}{b}{fn.loci}{fn.up}"
                  f"{fn.plot}", rt)
        safe_copy("fort.61",
                  f"{e_}{fn.dos}{b}{fn.tot}{fn.plot}.61", rt)
        safe_move("fort.61",
                  f"{ph}/{e_}{fn.dos}{b}{fn.tot}{fn.dn}"
                  f"{fn.plot}", rt)
        safe_copy("fort.71",
                  f"{e_}{fn.dos}{b}{fn.par}{fn.raw}.71", rt)
        safe_move("fort.71",
                  f"{ph}/{e_}{fn.dos}{b}{fn.par}{fn.dn}"
                  f"{fn.raw}", rt)
        safe_copy("fort.81",
                  f"{e_}{fn.dos}{b}{fn.loci}{fn.plot}.81",
                  rt)
        safe_move("fort.81",
                  f"{ph}/{e_}{fn.dos}{b}{fn.loci}{fn.dn}"
                  f"{fn.plot}", rt)
    else:
        safe_copy("fort.60",
                  f"{e_}{fn.dos}{b}{fn.tot}{fn.plot}.60", rt)
        safe_move("fort.60",
                  f"{ph}/{e_}{fn.dos}{b}{fn.tot}"
                  f"{fn.plot}", rt)
        safe_copy("fort.70",
                  f"{e_}{fn.dos}{b}{fn.par}{fn.raw}.70", rt)
        safe_move("fort.70",
                  f"{ph}/{e_}{fn.dos}{b}{fn.par}"
                  f"{fn.raw}", rt)
        safe_copy("fort.80",
                  f"{e_}{fn.dos}{b}{fn.loci}{fn.plot}.80",
                  rt)
        safe_move("fort.80",
                  f"{ph}/{e_}{fn.dos}{b}{fn.loci}"
                  f"{fn.plot}", rt)


def _manage_bond_output(e_, b, fn, ph, sp, rt):
    """Bond order + Q* output files."""
    if sp == 2:
        safe_copy("fort.10",
                  f"{e_}{fn.bond}{b}{fn.raw}.10", rt)
        safe_move("fort.10",
                  f"{ph}/{e_}{fn.bond}{b}{fn.up}"
                  f"{fn.raw}", rt)
        safe_copy("fort.11",
                  f"{e_}{fn.bond}{b}{fn.raw}.11", rt)
        safe_move("fort.11",
                  f"{ph}/{e_}{fn.bond}{b}{fn.dn}"
                  f"{fn.raw}", rt)
        if os.path.exists("fort.12"):
            safe_copy("fort.12",
                      f"{e_}{fn.bond}{b}{fn.bo3c}"
                      f"{fn.raw}.12", rt)
            safe_move("fort.12",
                      f"{ph}/{e_}{fn.bond}{b}{fn.bo3c}"
                      f"{fn.up}{fn.raw}", rt)
            safe_copy("fort.13",
                      f"{e_}{fn.bond}{b}{fn.bo3c}"
                      f"{fn.raw}.13", rt)
            safe_move("fort.13",
                      f"{ph}/{e_}{fn.bond}{b}{fn.bo3c}"
                      f"{fn.dn}{fn.raw}", rt)
    else:
        safe_copy("fort.10",
                  f"{e_}{fn.bond}{b}{fn.raw}.10", rt)
        safe_move("fort.10",
                  f"{ph}/{e_}{fn.bond}{b}{fn.raw}", rt)
        if os.path.exists("fort.12"):
            safe_copy("fort.12",
                      f"{e_}{fn.bond}{b}{fn.bo3c}"
                      f"{fn.raw}.12", rt)
            safe_move("fort.12",
                      f"{ph}/{e_}{fn.bond}{b}{fn.bo3c}"
                      f"{fn.raw}", rt)


def _manage_dimo_output(e_, b, fn, ph, sp, rt):
    """Dipole moment output files."""
    if sp == 2:
        safe_copy("fort.74",
                  f"{e_}{fn.dimo}{b}{fn.tot}"
                  f"{fn.plot}.74", rt)
        safe_move("fort.74",
                  f"{ph}/{e_}{fn.dimo}{b}{fn.tot}{fn.up}"
                  f"{fn.plot}", rt)
        safe_copy("fort.75",
                  f"{e_}{fn.dimo}{b}{fn.tot}"
                  f"{fn.plot}.61", rt)
        safe_move("fort.75",
                  f"{ph}/{e_}{fn.dimo}{b}{fn.tot}{fn.dn}"
                  f"{fn.plot}", rt)
    else:
        safe_copy("fort.74",
                  f"{e_}{fn.dimo}{b}{fn.tot}"
                  f"{fn.plot}.74", rt)
        safe_move("fort.74",
                  f"{ph}/{e_}{fn.dimo}{b}{fn.tot}"
                  f"{fn.plot}", rt)


def _manage_optc_output(e_, b, fn, ph, sp, rt):
    """Valence band optical properties output files."""
    if sp == 2:  # Spin polarized case.
        # In all cases, the total optc output file set from
        #   OLCAO_OPTC and OLCAOkkc needs to be copied.
        safe_copy("fort.40",
                  f"{e_}{fn.optc}{b}{fn.tot}{fn.cond}"
                  f"{fn.up}.40", rt)
        safe_copy("fort.50",
                  f"{e_}{fn.optc}{b}{fn.tot}{fn.eps2}"
                  f"{fn.up}.50", rt)
        safe_move("fort.40",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}{fn.cond}"
                  f"{fn.up}{fn.plot}", rt)
        safe_move("fort.50",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}{fn.eps2}"
                  f"{fn.up}{fn.plot}", rt)
        safe_move("fort.100",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}{fn.up}"
                  f"{fn.plot}", rt)
        safe_move("fort.110",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}{fn.eps1}"
                  f"{fn.up}{fn.plot}", rt)
        safe_move("fort.120",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}{fn.elf}"
                  f"{fn.up}{fn.plot}", rt)
        safe_move("fort.130",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.refrac}{fn.up}{fn.plot}", rt)
        safe_move("fort.140",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.extnct}{fn.up}{fn.plot}", rt)
        safe_move("fort.150",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.eps1i}{fn.up}{fn.plot}", rt)
        safe_move("fort.160",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.reflec}{fn.up}{fn.plot}", rt)
        safe_move("fort.170",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.absorp}{fn.up}{fn.plot}", rt)

        safe_copy("fort.41",
                  f"{e_}{fn.optc}{b}{fn.tot}{fn.cond}"
                  f"{fn.dn}.41", rt)
        safe_copy("fort.51",
                  f"{e_}{fn.optc}{b}{fn.tot}{fn.eps2}"
                  f"{fn.dn}.51", rt)
        safe_move("fort.41",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}{fn.cond}"
                  f"{fn.dn}{fn.plot}", rt)
        safe_move("fort.51",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}{fn.eps2}"
                  f"{fn.dn}{fn.plot}", rt)
        safe_move("fort.101",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}{fn.dn}"
                  f"{fn.plot}", rt)
        safe_move("fort.111",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}{fn.eps1}"
                  f"{fn.dn}{fn.plot}", rt)
        safe_move("fort.121",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}{fn.elf}"
                  f"{fn.dn}{fn.plot}", rt)
        safe_move("fort.131",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.refrac}{fn.dn}{fn.plot}", rt)
        safe_move("fort.141",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.extnct}{fn.dn}{fn.plot}", rt)
        safe_move("fort.151",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.eps1i}{fn.dn}{fn.plot}", rt)
        safe_move("fort.161",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.reflec}{fn.dn}{fn.plot}", rt)
        safe_move("fort.171",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.absorp}{fn.dn}{fn.plot}", rt)

        # Check for the existence of a certain file (fort.240)
        #   to decide if this optical properties calculation
        #   contains POPTC results.
        if os.path.exists("fort.240"):  # POPTC Scenario
            _manage_poptc_spin(e_, b, fn, ph, rt)
    else:  # Spin non-polarized case.
        # Always do the total optc results.
        safe_copy("fort.40",
                  f"{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.cond}.40", rt)
        safe_copy("fort.50",
                  f"{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.eps2}.50", rt)
        safe_move("fort.40",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}{fn.cond}"
                  f"{fn.plot}", rt)
        safe_move("fort.50",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}{fn.eps2}"
                  f"{fn.plot}", rt)
        safe_move("fort.100",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.plot}", rt)
        safe_move("fort.110",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}{fn.eps1}"
                  f"{fn.plot}", rt)
        safe_move("fort.120",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}{fn.elf}"
                  f"{fn.plot}", rt)
        safe_move("fort.130",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.refrac}{fn.plot}", rt)
        safe_move("fort.140",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.extnct}{fn.plot}", rt)
        safe_move("fort.150",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.eps1i}{fn.plot}", rt)
        safe_move("fort.160",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.reflec}{fn.plot}", rt)
        safe_move("fort.170",
                  f"{ph}/{e_}{fn.optc}{b}{fn.tot}"
                  f"{fn.absorp}{fn.plot}", rt)

        if os.path.exists("fort.240"):  # Partial optc
            _manage_poptc_nonspin(e_, b, fn, ph, rt)


def _manage_poptc_spin(e_, b, fn, ph, rt):
    """Partial optical properties output - spin polarized."""
    safe_copy("fort.209",
              f"{e_}{fn.optc}{b}{fn.par}{fn.kkc}"
              f"{fn.dat}.209", rt)
    safe_move("fort.209",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.kkc}"
              f"{fn.dat}", rt)

    # Spin up partial results.
    safe_copy("fort.240",
              f"{e_}{fn.optc}{b}{fn.par}{fn.cond}{fn.up}"
              f"{fn.raw}.240", rt)
    safe_move("fort.240",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.cond}"
              f"{fn.up}{fn.raw}", rt)
    safe_copy("fort.250",
              f"{e_}{fn.optc}{b}{fn.par}{fn.eps2}{fn.up}"
              f"{fn.raw}.250", rt)
    safe_move("fort.250",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.eps2}"
              f"{fn.up}{fn.raw}", rt)
    safe_move("fort.300",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.up}"
              f"{fn.raw}", rt)
    safe_move("fort.310",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.eps1}"
              f"{fn.up}{fn.raw}", rt)
    safe_move("fort.320",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.elf}"
              f"{fn.up}{fn.raw}", rt)
    safe_move("fort.330",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.refrac}"
              f"{fn.up}{fn.raw}", rt)
    safe_move("fort.340",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.extnct}"
              f"{fn.up}{fn.raw}", rt)
    safe_move("fort.350",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.eps1i}"
              f"{fn.up}{fn.raw}", rt)
    safe_move("fort.360",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.reflec}"
              f"{fn.up}{fn.raw}", rt)
    safe_move("fort.370",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.absorp}"
              f"{fn.up}{fn.raw}", rt)

    # Spin down partial results.
    safe_copy("fort.241",
              f"{e_}{fn.optc}{b}{fn.par}{fn.cond}{fn.dn}"
              f"{fn.raw}.241", rt)
    safe_move("fort.241",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.cond}"
              f"{fn.dn}{fn.raw}", rt)
    safe_copy("fort.251",
              f"{e_}{fn.optc}{b}{fn.par}{fn.eps2}{fn.dn}"
              f"{fn.raw}.251", rt)
    safe_move("fort.251",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.eps2}"
              f"{fn.dn}{fn.raw}", rt)
    safe_move("fort.301",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.dn}"
              f"{fn.raw}", rt)
    safe_move("fort.311",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.eps1}"
              f"{fn.dn}{fn.raw}", rt)
    safe_move("fort.321",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.elf}"
              f"{fn.dn}{fn.raw}", rt)
    safe_move("fort.331",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.refrac}"
              f"{fn.dn}{fn.raw}", rt)
    safe_move("fort.341",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.extnct}"
              f"{fn.dn}{fn.raw}", rt)
    safe_move("fort.351",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.eps1i}"
              f"{fn.dn}{fn.raw}", rt)
    safe_move("fort.361",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.reflec}"
              f"{fn.dn}{fn.raw}", rt)
    safe_move("fort.371",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.absorp}"
              f"{fn.dn}{fn.raw}", rt)


def _manage_poptc_nonspin(e_, b, fn, ph, rt):
    """Partial optical properties output - non-spin-pol."""
    safe_copy("fort.209",
              f"{e_}{fn.optc}{b}{fn.par}{fn.kkc}"
              f"{fn.dat}.209", rt)
    safe_move("fort.209",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.kkc}"
              f"{fn.dat}", rt)
    safe_copy("fort.240",
              f"{e_}{fn.optc}{b}{fn.par}{fn.cond}"
              f"{fn.raw}.240", rt)
    safe_move("fort.240",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.cond}"
              f"{fn.raw}", rt)
    safe_copy("fort.250",
              f"{e_}{fn.optc}{b}{fn.par}{fn.eps2}"
              f"{fn.raw}.250", rt)
    safe_move("fort.250",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.eps2}"
              f"{fn.raw}", rt)
    safe_move("fort.300",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}"
              f"{fn.raw}", rt)
    safe_move("fort.310",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.eps1}"
              f"{fn.raw}", rt)
    safe_move("fort.320",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.elf}"
              f"{fn.raw}", rt)
    safe_move("fort.330",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.refrac}"
              f"{fn.raw}", rt)
    safe_move("fort.340",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.extnct}"
              f"{fn.raw}", rt)
    safe_move("fort.350",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.eps1i}"
              f"{fn.raw}", rt)
    safe_move("fort.360",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.reflec}"
              f"{fn.raw}", rt)
    safe_move("fort.370",
              f"{ph}/{e_}{fn.optc}{b}{fn.par}{fn.absorp}"
              f"{fn.raw}", rt)


def _manage_pacs_output(e_, b, fn, ph, sp, rt):
    """Photo-absorption cross section output files."""
    if sp == 2:
        safe_copy("fort.50",
                  f"{e_}{fn.pacs}{b}{fn.plot}{fn.up}.50",
                  rt)
        safe_move("fort.50",
                  f"{ph}/{e_}{fn.pacs}{b}{fn.up}"
                  f"{fn.plot}", rt)
        safe_copy("fort.51",
                  f"{e_}{fn.pacs}{b}{fn.plot}{fn.dn}.51",
                  rt)
        safe_move("fort.51",
                  f"{ph}/{e_}{fn.pacs}{b}{fn.dn}"
                  f"{fn.plot}", rt)
    else:
        safe_copy("fort.50",
                  f"{e_}{fn.pacs}{b}{fn.plot}.50", rt)
        safe_move("fort.50",
                  f"{ph}/{e_}{fn.pacs}{b}{fn.plot}", rt)


def _manage_nlop_output(e_, b, fn, ph, sp, rt):
    """Non-linear optical properties output files."""
    if sp == 2:
        safe_copy("fort.50",
                  f"{e_}{fn.optc}{b}{fn.chi2}{fn.up}.50",
                  rt)
        safe_move("fort.50",
                  f"{ph}/{e_}{fn.optc}{b}{fn.chi2}{fn.up}"
                  f"{fn.plot}", rt)
        safe_move("fort.100",
                  f"{ph}/{e_}{fn.optc}{b}{fn.up}"
                  f"{fn.plot}", rt)
        safe_move("fort.110",
                  f"{ph}/{e_}{fn.optc}{b}{fn.chi1}{fn.up}"
                  f"{fn.plot}", rt)

        safe_copy("fort.51",
                  f"{e_}{fn.optc}{b}{fn.chi2}{fn.dn}.51",
                  rt)
        safe_move("fort.51",
                  f"{ph}/{e_}{fn.optc}{b}{fn.chi2}{fn.dn}"
                  f"{fn.plot}", rt)
        safe_move("fort.101",
                  f"{ph}/{e_}{fn.optc}{b}{fn.dn}"
                  f"{fn.plot}", rt)
        safe_move("fort.111",
                  f"{ph}/{e_}{fn.optc}{b}{fn.chi1}{fn.dn}"
                  f"{fn.plot}", rt)
    else:
        safe_copy("fort.50",
                  f"{e_}{fn.optc}{b}{fn.chi2}.50", rt)
        safe_move("fort.50",
                  f"{ph}/{e_}{fn.optc}{b}{fn.chi2}"
                  f"{fn.plot}", rt)
        safe_move("fort.100",
                  f"{ph}/{e_}{fn.optc}{b}{fn.plot}", rt)
        safe_move("fort.110",
                  f"{ph}/{e_}{fn.optc}{b}{fn.chi1}"
                  f"{fn.plot}", rt)


def _manage_sige_output(e_, b, fn, ph, sp, rt):
    """Sigma(E) output files."""
    if sp == 2:
        safe_copy("fort.50",
                  f"{e_}{fn.sige}{b}{fn.cond}{fn.up}.50",
                  rt)
        safe_move("fort.50",
                  f"{ph}/{e_}{fn.sige}{b}{fn.cond}{fn.up}"
                  f"{fn.plot}", rt)
        safe_copy("fort.51",
                  f"{e_}{fn.sige}{b}{fn.cond}{fn.dn}.51",
                  rt)
        safe_move("fort.51",
                  f"{ph}/{e_}{fn.sige}{b}{fn.cond}{fn.dn}"
                  f"{fn.plot}", rt)
    else:
        safe_copy("fort.50",
                  f"{e_}{fn.sige}{b}{fn.cond}.50", rt)
        safe_move("fort.50",
                  f"{ph}/{e_}{fn.sige}{b}{fn.cond}"
                  f"{fn.plot}", rt)


def _manage_sybd_output(e_, b, fn, ph, sp, rt):
    """Symmetric band structure output files."""
    if sp == 2:
        safe_move("fort.31",
                  f"{e_}{fn.sybd}{b}{fn.raw}.31", rt)
        safe_copy("fort.33",
                  f"{ph}/{fn.vdim}{b}{fn.raw}", rt)
        safe_move("fort.33",
                  f"{e_}{fn.sybd}{b}{fn.raw}.33", rt)
        safe_move("fort.41",
                  f"{ph}/{e_}{fn.sybd}{b}{fn.up}"
                  f"{fn.plot}", rt)
        safe_move("fort.32",
                  f"{e_}{fn.sybd}{b}{fn.raw}.32", rt)
        safe_copy("fort.34",
                  f"{ph}/{fn.vdim}{b}{fn.raw}", rt)
        safe_move("fort.34",
                  f"{e_}{fn.sybd}{b}{fn.raw}.34", rt)
        safe_move("fort.42",
                  f"{ph}/{e_}{fn.sybd}{b}{fn.dn}"
                  f"{fn.plot}", rt)
    else:
        safe_move("fort.31",
                  f"{e_}{fn.sybd}{b}{fn.raw}.31", rt)
        safe_copy("fort.33",
                  f"{ph}/{fn.vdim}{b}{fn.raw}", rt)
        safe_move("fort.33",
                  f"{e_}{fn.sybd}{b}{fn.raw}.33", rt)
        safe_move("fort.41",
                  f"{ph}/{e_}{fn.sybd}{b}{fn.plot}", rt)


def _manage_force_output(e_, b, fn, ph, sp, rt):
    """Force output files."""
    if sp == 2:
        safe_copy("fort.98",
                  f"{e_}{fn.force}{b}{fn.dat}.98", rt)
        safe_move("fort.98",
                  f"{ph}/{e_}{fn.force}{b}{fn.up}"
                  f"{fn.dat}", rt)
        safe_copy("fort.99",
                  f"{e_}{fn.force}{b}{fn.dat}.99", rt)
        safe_move("fort.99",
                  f"{ph}/{e_}{fn.force}{b}{fn.dn}"
                  f"{fn.dat}", rt)
    else:
        safe_copy("fort.98",
                  f"{e_}{fn.force}{b}{fn.dat}.98", rt)
        safe_move("fort.98",
                  f"{ph}/{e_}{fn.force}{b}{fn.dat}", rt)


def _manage_field_output(e_, b, fn, ph, sp, rt):
    """Field output files (wave function, charge density,
    potential, profiles, charge centers)."""
    # If the profile files were created, then move them.
    if os.path.exists("fort.30"):
        safe_move("fort.30",
                  f"{ph}/{e_}{fn.field}{b}{fn.profile}"
                  f"-a{fn.dat}", rt)
        safe_move("fort.31",
                  f"{ph}/{e_}{fn.field}{b}{fn.profile}"
                  f"-b{fn.dat}", rt)
        safe_move("fort.32",
                  f"{ph}/{e_}{fn.field}{b}{fn.profile}"
                  f"-c{fn.dat}", rt)

    # If the charge center file was created, then move it.
    if os.path.exists("fort.56"):
        safe_move("fort.56",
                  f"{ph}/{e_}{fn.field}{b}{fn.rho}"
                  f"{fn.dat}", rt)

    # Only move XDMF-HDF5 files if they were created.
    if os.path.exists("fort.78"):
        safe_move("fort.78",
                  f"{ph}/{e_}{fn.field}{b}{fn.xdmf3}", rt)


def _manage_mtop_output(e_, b, fn, ph, sp, rt):
    """Polarization (mtop) output files."""
    if sp == 2:
        safe_copy("fort.180",
                  f"{e_}{fn.mtop}{b}{fn.tot}"
                  f"{fn.plot}.180", rt)
        safe_move("fort.180",
                  f"{ph}/{e_}{fn.mtop}{b}{fn.tot}{fn.up}"
                  f"{fn.plot}", rt)
        safe_copy("fort.181",
                  f"{e_}{fn.mtop}{b}{fn.tot}"
                  f"{fn.plot}.181", rt)
        safe_move("fort.181",
                  f"{ph}/{e_}{fn.mtop}{b}{fn.tot}{fn.dn}"
                  f"{fn.plot}", rt)
    else:
        safe_copy("fort.180",
                  f"{e_}{fn.mtop}{b}{fn.tot}"
                  f"{fn.plot}.180", rt)
        safe_move("fort.180",
                  f"{ph}/{e_}{fn.mtop}{b}{fn.tot}"
                  f"{fn.plot}", rt)


def _manage_loen_output(e_, b, fn, ph, rt):
    """Local environment output files."""
    # Move the local environment structural analysis files
    #   if created.
    if os.path.exists("fort.21"):
        safe_copy("fort.21",
                  f"{e_}{fn.loen}{b}{fn.plot}.21", rt)
        safe_move("fort.21",
                  f"{ph}/{e_}{fn.loen}{b}{fn.plot}", rt)


# ------------------------------------------------------------------ #
#                       PACS Update                                   #
# ------------------------------------------------------------------ #

# Hartree energy in eV (from structure_control.py constants).
HARTREE_EV = 27.21138386


def update_pacs(fin_scf_edge, update_scf_basis, fn,
                proj_home, runtime_fh):
    """Update the PACS input data in olcao.dat after both SCF
    calculations (initial and final state) are complete.

    This subroutine determines the approximate location of the
    edge onset based on the total energy difference between
    the two states (initial/final). It then replaces any
    "UNKNOWN" (-1) values in the PACS_INPUT_DATA section of
    olcao.dat with the computed total energy difference.

    Args:
        fin_scf_edge: The final state edge string (e.g. "1s").
        update_scf_basis: The SCF basis string (e.g. "fb").
        fn: FileNames instance.
        proj_home: Project home directory.
        runtime_fh: Open RUNTIME file handle.
    """
    edge_ = f"{fin_scf_edge}_"
    basis = f"-{update_scf_basis}"
    dat = fn.dat
    out = fn.out
    scf = fn.scf
    olcao = fn.olcao

    orig_dir = os.getcwd()
    os.chdir(proj_home)

    # Read the initial state SCF output to find total energy.
    init_scf_path = f"gs_{scf}{basis}{out}"
    if not os.path.exists(init_scf_path):
        sys.exit(
            f"Cannot open {init_scf_path} for reading."
        )

    init_energy = None
    with open(init_scf_path, 'r') as f:
        for line in f:
            if "TOTAL ENERGY" in line:
                parts = line.split()
                init_energy = float(parts[-1])

    # Read the final state SCF output to find total energy.
    fin_scf_path = f"{edge_}{scf}{basis}{out}"
    if not os.path.exists(fin_scf_path):
        sys.exit(
            f"Cannot open {fin_scf_path} for reading."
        )

    fin_energy = None
    with open(fin_scf_path, 'r') as f:
        for line in f:
            if "TOTAL ENERGY" in line:
                parts = line.split()
                fin_energy = float(parts[-1])

    # Compute the total energy difference (TEDiff).
    te_diff = abs(fin_energy - init_energy) * HARTREE_EV

    # Read the olcao.dat file, replace UNKNOWN values, and
    #   write to a temp file.
    olcao_dat = f"{olcao}{dat}"
    if not os.path.exists(olcao_dat):
        sys.exit(
            f"Cannot open {olcao_dat} for reading."
        )

    unknown_found = False
    with open(olcao_dat, 'r') as f_in, \
            open("pacs_temp", 'w') as f_out:
        for line in f_in:
            if "PACS_INPUT_DATA" in line:
                f_out.write(line)
                # Write next 4 lines as-is.
                for _ in range(4):
                    next_line = next(f_in)
                    f_out.write(next_line)
                # The 4th line has the number of core
                #   orbitals.
                parts = next_line.split()
                num_core_orbitals = int(parts[0])

                for _ in range(num_core_orbitals):
                    orb_line = next(f_in)
                    values = orb_line.split()
                    if (float(values[4]) == -1
                            and int(values[0])
                            == EDGE_MAP.get(
                                fin_scf_edge, (0, 0)
                            )[0]
                            and int(values[1])
                            == EDGE_MAP.get(
                                fin_scf_edge, (0, 0)
                            )[1]):
                        f_out.write(
                            f"{values[0]} {values[1]} "
                            f"{values[2]} {values[3]} "
                            f"{te_diff} {values[5]}"
                            f"  # QN_n QN_l Init1 Init2"
                            f" TEDiff\n"
                        )
                        unknown_found = True
                    else:
                        f_out.write(orb_line)
            else:
                f_out.write(line)

    # Update the version of the pacs input file.
    safe_move("pacs_temp", olcao_dat, runtime_fh)

    # Only print a notice if we changed anything.
    if unknown_found:
        runtime_fh.write(f"{olcao_dat} Updated.\n")

    os.chdir(orig_dir)


# ------------------------------------------------------------------ #
#                          Cleanup                                    #
# ------------------------------------------------------------------ #

def clean_up(temp, fn, runtime_fh):
    """Clean up the lock file and close the runtime log."""
    lock_path = os.path.join(temp, fn.olcao_lock)
    if os.path.exists(lock_path):
        os.remove(lock_path)
    runtime_fh.write("Program Sequence Complete.\n")
    runtime_fh.close()


# ------------------------------------------------------------------ #
#                         Main Program                                #
# ------------------------------------------------------------------ #

def main():

    # Get script settings from a combination of the resource
    #   control file and parameters given by the user on the
    #   command line.
    settings = ScriptSettings()

    # Initialize file name components.
    fn = FileNames()

    # Set up the execution environment.
    # Initialize the directories to be used.
    proj_home, inputs, temp = init_directories(fn)

    # Initialize the IO for storing results to RUNTIME.
    runtime_fh = init_io(proj_home, fn)

    # Initialize the executables and execution method.
    executable, exe_mechanism, sub_mechanism = init_exes(
        settings, fn, inputs
    )
    bin_dir = os.environ.get('OLCAO_BIN', '')

    # Determine the XC Code and if the calculation is spin
    #   polarized. Use 1 for non-spinpol; use 2 for spinpol.
    spin_pol = xc_code_spin(proj_home, inputs, fn)

    # Check that a lock file does not exist blocking
    #   execution. If it doesn't, then make one so we can
    #   start running the job.
    lock_path = os.path.join(temp, fn.olcao_lock)
    if os.path.exists(lock_path):
        olcao_exit(
            f"Lock file found in {temp}.\n"
            f"Is another olcao script running?\n"
            f"Did an olcao script die badly?\n",
            runtime_fh,
        )
    else:
        # Create the lock file.
        with open(lock_path, 'w') as f:
            pass

    # --- Run the requested job --- #

    # Change the working directory to the intermediate
    #   location.
    os.chdir(temp)

    # Mark the date and time of the beginning of this run.
    output_label = (
        f"{settings.edge}_{settings.job_name}"
        f"-{settings.basis_scf}-{settings.basis_pscf}"
    )
    runtime_fh.write(
        f"{output_label}----spin={spin_pol}"
        f"-----------------------------\n"
    )
    runtime_fh.write(
        f"Start:  {datetime.now()}\n"
    )

    # Mark the lock file with the current executable.
    with open(lock_path, 'w') as f:
        f.write(f"{output_label}\n")

    # Check for and stage the necessary input files.
    manage_input(
        settings, fn, proj_home, inputs, temp, runtime_fh
    )

    # Create the proper command line parameters (CLP) for
    #   this job.
    serial_flag = 1 if settings.serialxyz else 0
    job_clp = (
        f"{settings.basis_code_scf} "
        f"{settings.basis_code_pscf} "
        f"{settings.qn_n} {settings.qn_l} "
        f"{settings.job_id} {serial_flag}"
    )

    # Execute the requested task program and save the
    #   runtime data.
    execute_program(
        job_clp, settings, fn, bin_dir,
        executable, exe_mechanism, sub_mechanism,
        spin_pol, runtime_fh,
    )

    # Manage the output files.
    manage_output(
        settings, fn, proj_home, spin_pol, runtime_fh
    )

    # Mark the date and time of the ending of this task.
    runtime_fh.write(
        f"End:  {datetime.now()}\n"
    )

    # Check for a kill file request.
    if os.path.exists(fn.olcao_kill):
        runtime_fh.close()
        sys.exit()

    # Clean up any leftovers and tie up loose ends.
    clean_up(temp, fn, runtime_fh)


if __name__ == '__main__':
    # Everything before this point was a subroutine definition
    #   or a request to import information from external
    #   modules. Only now do we actually start running the
    #   program. The purpose of this is to allow another
    #   python program to import *this* script and call its
    #   functions internally.
    main()
