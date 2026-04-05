#!/usr/bin/env python3

"""makeinput.py -- Generate all OLCAO input files from an olcao.skl skeleton.

This is the Python port of the Perl ``makeinput`` script.  It reads an
olcao.skl file and produces the full set of input files needed for any
OLCAO calculation (SCF, band, DOS, bond, optical, XANES, etc.).

The script follows the XYZ.py / XYZrc.py pattern:
  1. Load defaults from makeinputrc.py ($OLCAO_RC or cwd).
  2. Parse the command line with argparse.
  3. Reconcile CLI args with rc defaults.
  4. Execute the main workflow.


CONSIDERATIONS AS YOU PREPARE YOUR INPUT FILE
==============================================

(1) Please remember the following definitions:
    a. Elements are defined by the periodic table of the elements.
    b. Species are defined by the structure of the system (local or global).
    c. Types are defined by the needs of the calculation.

(2) One of the key components of running an efficient and accurate
    calculation is to properly assign the potential types of the atoms in
    the system.  In many (most?) cases the potential types and the species
    are identical sets of atoms.  When they differ it is usually done
    automatically by the makeinput script for specific calculations (e.g.
    ELNES/XANES).

(3) You will assign the types depending on the kind of system you are trying
    to run a calculation on.
    a. Is the system crystalline?  If so, types will likely be assigned by
       crystallographic equivalency within the olcao.skl input file.  This
       makeinput script will then create atoms in the appropriate positions
       based on the space group and equivalent atoms will have the same
       potential type.  There is little to consider.
    b. Is the system crystalline with a small perturbation?  (e.g. atom
       substitution, vacancy, interstitial, etc.)  This can be treated with
       a two step process.  In the first step the crystalline olcao.skl file
       is created and the makeinput script is called.  This will produce a
       set of input files in the ./inputs directory.  You will then edit the
       file: ./inputs/olcao.mi (makeinput) which contains all atoms
       explicitly listed.  You can apply your perturbation(s) and then copy
       the ./inputs/olcao.mi file back out as the new input file.  The
       species can be assigned using the 'reduce', 'target', or 'block'
       methods as the perturbation nature may require.
    c. Is the system a complex system with no symmetry or is it a crystalline
       system after a relaxation (done by some other program)?  You will
       likely need to use the 'reduce' scheme to assign types.

(4) What follows is a short review of the different grouping techniques.
    These assignment methods can be overlapped so that a system can be
    defined by the olcao.skl file and then further modified by xanes
    requests or other targeted requests or similarity classifications etc.
    a. Atoms are initially explicitly grouped in the olcao.skl file.
    b. Atoms can be grouped by similarity characteristics that detail the
       threshold where two atoms are considered to be in the same group
       based on their local configuration.  (Reduce scheme)
    c. Atoms can be grouped by a list of targeted locations and radius
       cutoffs.  Then, either inside or outside the cutoff radius, atoms
       can be grouped alike or unalike.  (Target scheme)
    d. Atoms can be grouped by blocks where inside or outside the specified
       3D lattice-like boundaries atoms will be grouped alike or unalike.
       (Block scheme)
    e. Atoms can be grouped in a special way that is related to the targeted
       method given above.  For XANES calculations the targets can be either
       explicitly requested or automatically selected from each species.
       For each target a separate input file set is created so the targets
       do not overlap.  Within a sphere around the target all atoms will be
       given different types and the target atom will have its core orbitals
       included in the calculation instead of being orthogonalized out.
       This is just an automation of the target method for easily creating
       lots of input files for a given system for XAS calculation.
       (XANES Target scheme)
    f. Whatever future ideas I can come up with that are needed and useful.
"""

import argparse as ap
import math
import os
import shutil
import sys
from datetime import datetime


# ---------------------------------------------------------------------------
# ScriptSettings -- command-line parameters and rc-file defaults
# ---------------------------------------------------------------------------

class ScriptSettings:
    """Holds all user-controllable parameters for makeinput.

    On construction the object:
      1. Reads defaults from makeinputrc.py.
      2. Parses the command line (argparse).
      3. Reconciles CLI values with rc defaults.
      4. Records the command line to the ``command`` file.
    """

    def __init__(self):
        # Read default variables from the resource control file.
        rc = self._load_rc()
        self.assign_rc_defaults(rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the rc defaults.
        self.reconcile(args)

        # Record the command line that was used.
        self.record_clp()

    # ------------------------------------------------------------------
    # RC loading
    # ------------------------------------------------------------------

    @staticmethod
    def _load_rc():
        """Import makeinputrc and return its parameter dictionary.

        Search order: cwd, then $OLCAO_RC.
        """
        # Try cwd first for a local override.
        cwd_rc = os.path.join(os.getcwd(), "makeinputrc.py")
        if os.path.isfile(cwd_rc):
            sys.path.insert(0, os.getcwd())
        else:
            rc_dir = os.getenv("OLCAO_RC")
            if not rc_dir:
                sys.exit("Error: $OLCAO_RC is not set and no local "
                         "makeinputrc.py found. See instructions.")
            sys.path.insert(0, rc_dir)

        from makeinputrc import parameters_and_defaults
        return parameters_and_defaults()

    # ------------------------------------------------------------------
    # Assign rc defaults
    # ------------------------------------------------------------------

    def assign_rc_defaults(self, rc):
        """Populate instance attributes from the rc dictionary."""

        # Executable names.
        self.makekpoints_exec = rc["makekpoints_exec"]
        self.contract_exec = rc["contract_exec"]

        # Potential modification.
        self.mod_pot = rc["mod_pot"]
        self.mod_element_name = rc["mod_element_name"]
        self.min_mod_term = rc["min_mod_term"]
        self.max_mod_term = rc["max_mod_term"]
        self.num_mod_terms = rc["num_mod_terms"]

        # Output file flags.
        self.pdb = rc["pdb"]
        self.cif = rc["cif"]
        self.queue_type = rc["queue_type"]

        # Cutoffs.
        self.bf_cutoff = rc["bf_cutoff"]
        self.es_cutoff = rc["es_cutoff"]

        # Exchange correlation mesh.
        self.xc_code = rc["xc_code"]
        self.xc_in_weight = rc["xc_in_weight"]
        self.xc_out_weight = rc["xc_out_weight"]
        self.xc_in_samp = rc["xc_in_samp"]
        self.xc_out_samp = rc["xc_out_samp"]
        self.xc_spacing_samp = rc["xc_spacing_samp"]
        self.num_samp_vectors = rc["num_samp_vectors"]

        # DOS.
        self.e_delta_dos = rc["e_delta_dos"]
        self.sigma_dos = rc["sigma_dos"]
        self.emin_dos = rc["emin_dos"]
        self.emax_dos = rc["emax_dos"]
        self.detail_code_pdos = rc["detail_code_pdos"]
        self.iter_tdos = rc["iter_tdos"]

        # BOND.
        self.max_len_bond = rc["max_len_bond"]
        self.emin_bond = rc["emin_bond"]
        self.emax_bond = rc["emax_bond"]
        self.sigma_bond = rc["sigma_bond"]
        self.e_delta_bond = rc["e_delta_bond"]
        self.output_bondq = rc["output_bondq"]
        self.bond_3c = rc["bond_3c"]
        self.max_neighbor_bond = rc["max_neighbor_bond"]

        # PACS.
        self.e_delta_pacs = rc["e_delta_pacs"]
        self.sigma_pacs = rc["sigma_pacs"]
        self.onset_slack_pacs = rc["onset_slack_pacs"]
        self.energy_window_pacs = rc["energy_window_pacs"]

        # OPTC.
        self.e_delta_optc = rc["e_delta_optc"]
        self.e_delta_nlop = rc["e_delta_nlop"]
        self.e_delta_sige = rc["e_delta_sige"]
        self.sigma_optc = rc["sigma_optc"]
        self.sigma_nlop = rc["sigma_nlop"]
        self.sigma_sige = rc["sigma_sige"]
        self.e_trans_optc = rc["e_trans_optc"]
        self.e_trans_nlop = rc["e_trans_nlop"]
        self.e_trans_sige = rc["e_trans_sige"]
        self.e_cutoff_optc = rc["e_cutoff_optc"]
        self.e_cutoff_nlop = rc["e_cutoff_nlop"]
        self.e_cutoff_sige = rc["e_cutoff_sige"]
        self.detail_code_poptc = rc["detail_code_poptc"]

        # MAIN (SCF).
        self.num_iter_main = rc["num_iter_main"]
        self.converg_main = rc["converg_main"]
        self.spin_main = rc["spin_main"]
        self.therm_smear_main = rc["therm_smear_main"]

        # SYBD / cell name.
        self.sybd_path = rc["sybd_path"]
        self.cell_name = rc["cell_name"]

        # K-points.
        self.kp_mesh_scf = list(rc["kp_mesh_scf"])
        self.kp_mesh_pscf = list(rc["kp_mesh_pscf"])
        self.kp_weight_sum = rc["kp_weight_sum"]
        self.kp_shift = rc["kp_shift"]
        self.print_bz = rc["print_bz"]
        self.scale_factor = rc["scale_factor"]
        self.kp_note = ["", "(General)", "(General)"]  # [0]=unused, [1]=scf, [2]=pscf
        # Density-based k-point mode.  When use_kp_density is
        # True, the makeKPoints program is bypassed and style-
        # code-2 k-point files are written directly.  The
        # density values are stored per group (1=SCF, 2=PSCF).
        self.use_kp_density = False
        self.kp_density = [0.0, 0.0, 0.0]  # [0]=unused

        # K-point integration method code per group.
        # 0 = Gaussian/histogram (default), 1 = LAT
        # (linear analytic tetrahedron).  Future methods
        # may use higher integers.  Index convention:
        # [0]=unused, [1]=SCF, [2]=PSCF.
        self.kp_intg_code = [0, 0, 0]

        # Basis / potential substitution lists (populated by CLI).
        self.basis_sub_out = []
        self.basis_sub_in = []
        self.pot_sub_out = []
        self.pot_sub_in = []

        # Grouping methods -- accumulated from CLI, each entry is a dict.
        self.methods = []       # ordered list of ("reduce"|"target"|"block", idx)
        self.reduces = []       # list of dicts, one per -reduce invocation
        self.targets = []       # list of dicts, one per -target invocation
        self.blocks = []        # list of dicts, one per -block invocation

        # XANES.
        self.xanes = rc["xanes"]
        self.xanes_radius = rc["xanes_radius"]
        self.xanes_atoms = []   # explicit atom numbers (empty = auto-select)

        # Misc.
        self.state_factor = rc["state_factor"]
        self.rel = rc["rel"]
        self.do_basis_vis = rc["do_basis_vis"]
        self.emu = rc["emu"]
        self.no_core = rc["no_core"]

        # Slurm.
        self.partition = rc["partition"]
        self.account = rc["account"]
        self.time = rc["time"]
        self.memory = rc["memory"]
        self.cpus = rc["cpus"]
        self.nodes = rc["nodes"]

    # ------------------------------------------------------------------
    # Command-line parsing
    # ------------------------------------------------------------------

    def parse_command_line(self):
        """Build the argparse parser, add all arguments, and parse."""

        prog_name = "makeinput"

        description_text = """\
Generate all OLCAO input files from an olcao.skl skeleton file.

Creates olcao.dat, structure.dat, scfV.dat, k-point files, basis/potential
data, and (optionally) XANES file sets, slurm scripts, PDB/CIF files, etc.

See the module docstring (top of makeinput.py) for important conceptual
background on elements, species, types, and the different grouping schemes
(reduce, target, block, XANES target).
"""

        epilog_text = """\
OPTIONS EXPLANATIONS
====================

-basisdb    Supply the location of an alternate atomic basis functions
            database.  The default value is $OLCAO_DATA/atomicBDB.

-potdb      Supply the location of an alternate potential function database.
            The default value is $OLCAO_DATA/atomicPDB.

-subbasis   When the input files are constructed, the value in SUBBASIS_IN
            will substitute for the value in SUBBASIS_OUT.  For example:
            "-subbasis b 2" will use the contract2.dat basis definition file
            instead of the default contract1.dat basis definition file for
            all boron.  This allows one to use alternate configurations from
            the default basis set.  The exact way the substitution is done
            depends on the value of SUBBASIS_OUT:  If it is an element (e.g.
            si, c, zr) then all atoms of that element will be substituted.
            If it is a species (e.g. si1, c3, zr2) then only that specific
            species will be substituted.  This substitution is for the basis
            set only.  This option is repeatable.

-subpot     Same behavior as -subbasis except that the potential is
            substituted instead of the basis.  This option is repeatable.

-modpot     Modify the potential for the element specified by ELEM so that
            the minimum value is MIN, the maximum value is MAX, and the
            number of terms is NTERMS.  Note that this option can only be
            applied to one element as the program is presently written.

-scfkp      Specify the mesh of k-points to be used for the SCF calculation.
            The three parameters are the size of the mesh in the a, b, and c
            directions.  (e.g. 2 2 2 would produce 8 kpoints before symmetry
            reductions, and 2 3 4 would produce 24 kpoints before symmetry
            reductions.)  This applies to the primary SCF programs (setup and
            main) and also to any subroutines run within the SCF stage (e.g.
            dos or bond).  Please note that most of the time things like dos
            or bond are automatically run in the post-SCF stage where a
            larger number of kpoints are typically used.  One further
            important note:  If the kpoint scheme given is 1 1 1 then it is
            assumed that the one point is at the gamma site.  If no kpoint
            definition is given then 1 general kpoint is used.  This is
            important since the program will run differently (faster) if the
            gamma kpoint is used since all the integral matrices will be real
            with no imaginary component.

-pscfkp     Same as -scfkp except that it applies to the post-SCF
            calculations (band, bond, dos, optc).  NOTE that sybd will
            choose kpoints given by the path associated with the crystalline
            cell or the -sybdpath parameter if it was given to request an
            explicit path.

-kp         Apply the same kpoint selection to the SCF and post-SCF kpoints.

-scfkpd     Specify the k-point volume density for the SCF
            calculation as a single number (kpoints per unit
            reciprocal-space volume, in Bohr^-3).  The total
            kpoint count will be density * V_BZ, distributed
            as uniformly as possible across the three axes.
            Instead of giving explicit mesh dimensions, the
            program writes a style-code-2 k-point file and
            lets OLCAO compute the per-axis mesh counts from
            the reciprocal cell geometry at runtime.  This
            bypasses the makeKPoints program entirely.

-pscfkpd    Same as -scfkpd but for the post-SCF calculations
            (band, bond, dos, optc).

-kpd        Apply the same k-point volume density to both SCF
            and post-SCF.  If any density option is given,
            density mode is used for both groups (an unset
            group defaults to density 1).  Density options
            are mutually exclusive with explicit mesh options
            (-kp, -scfkp, -pscfkp); if both are given,
            density mode wins and a warning is printed.

-scfkpint   Select the k-point integration method for the SCF
            calculation.  M is an integer: 0 = Gaussian/histogram
            (default), 1 = LAT (linear analytic tetrahedron).
            Higher integers are reserved for future methods.
            This value is written into the KPOINT_INTG_CODE field
            of the k-point input file.  Non-Gaussian methods
            currently require density mode (-kpd/-scfkpd/-pscfkpd);
            the legacy makeKPoints program does not support them.

-pscfkpint  Same as -scfkpint but for post-SCF calculations
            (band, bond, dos, optc).

-kpint      Apply the same integration method to both SCF and
            post-SCF.  Per-group flags (-scfkpint, -pscfkpint)
            override this if also given.

-kpshift    Force a specific shift to the kpoint mesh by the given fractional
            a, b, c amounts.

-printbz    Create an additional file as part of the input containing a
            description of the Brillouin zone suitable for importing in a
            Python script.  The two required options are an integer
            identifying which Brillouin zone to print and a scale factor
            that causes the size of the Brillouin zone to be scaled.

-xccode     Set the exchange correlation code (see $OLCAO_DIR/share/
            xc_code.dat for the list of codes).

-xcmesh     Define the construction of the real space mesh for sampling the
            charge density for the exchange correlation evaluation.  The mesh
            is spherical and atom centered.  Sub-options:
              -numvect    Number of radial sample vector directions.
              -weight     Weighting of points within/outside a radial cutoff.
              -samp       Distribution of points along the vector directions:
                          in-sampling, out-sampling, and spacing.

-target     Consider a point given either by an x,y,z location (-atxyz), an
            a,b,c location (-atabc), or an atom number (-atom).  Then the
            atoms either in or out of the zone (sphere of given -sphere
            radius) will be considered in terms of their -operand
            (species, type, reduce) and will be grouped by their -relate
            relation (alike or diff).  If the reduce operand is given then
            this target is used as a selecting tool for a reduce call.

-block      Functions just like the -target option except that the area is
            defined like a slab or block of the system using -abc to specify
            from/to values for each lattice direction.  Note that it is
            possible to use the letters 'a', 'b', or 'c' in the place of
            actual numbers when one wants the block "To" value to be the
            maximum.  (e.g. "-abc 0 a 0 b 0 c" would include the entire
            cell.)

-reduce     Collect information about each atom in order to group all the
            atoms based on the similarity of the information.  A group is
            defined by a series of spherical shells (levels), what neighbor
            atoms are in the shells, and what the distance is to each shell
            from the center.  Parameters can be adjusted to determine the
            number of shells (-level), the thickness of the shells (-thick),
            a maximum cutoff radius (-cutoff), and a distance threshold for
            comparing the same shell number from different atoms (-tolerance).
            The basic idea is to find the nearest atom to the current atom,
            define a shell of a given thickness and record what atoms are in
            that shell.  Then repeat (reduceLevel) times with the next
            nearest atom outside the shell.  Once all the shells have been
            defined for each atom, the results are compared in terms of the
            atoms in each shell and the shell distances.  This grouping
            method is not at all defined by boundaries and so does not have
            in and out -zone parameters.  It will also not have the ability
            to make grouping dissimilar (e.g. it will not find a group of
            atoms that are similar and then make them all different).  This
            method will only work to make species out of elements, and will
            not make types out of species.  It can be applied differently to
            different groups of atoms via -selection.

-xanes      Instead of generating one set of input files, make one set for
            each of the atomic species or listed xanes atoms in the system.
            Each input file set will be different in that a xanes atom of
            the current input file set will have the core orbitals included
            in the calculation, and the atoms within a given radius (default
            3.50 A) of the xanes atom will all have different types.
            Sub-options:
              -sphere     Radius around the target atom.
              -atom       List of target atom numbers (supports N..M ranges).

-sybdpath   Specify the type of cell and particular path to be used in the
            symmetric band calculation (SYBD).  You should be consistent with
            the actual cell type or you will get a warning (even though it
            will still allow you to create the input files).  The valid
            options are present in the $OLCAO_DATA/sybdDB directory.  Simply
            specify one of those names as the sybdPath.

-rel        Prepare the OLCAO input files for a relativistic calculation.
            (Not yet operational in the Fortran programs.)

-statefactor
            Select the number of states to calculate as a multiple of the
            number of valence electrons.  Default: 2.5.

-pdb        Generate a PDB (Protein Data Bank) crystal structure file.

-cif        Generate a CIF (Crystallographic Interchange Format) structure
            file.  This format allows for greater numbers of species for each
            atom and is the preferred method.

-basisVis   Ask the contract program to produce files in the .inputTemp
            directory for visualizing the complete atomic orbital basis.  The
            files are numerical data of the radial part of the basis
            functions, and a set of POVRay scene files (one for each type of
            atomic orbital in the system).

-nocore     Include core orbitals in the valence (no orthogonalization).
            This causes the basis functions that are normally orthogonalized
            out of the problem to instead be included as part of the list of
            valence basis functions so that they are also part of the SCF
            process.  This allows one to compute the core energy eigenvalues.
            If this option is not given, then basis functions that are
            designated as "core" will be orthogonalized out of the
            eigenvalue problem.

-emu        Generate EMU configuration files.

-slurm      Specify parameters for the slurm submission script.
            Sub-options:
              -p    Partition name.
              -a    Account to charge.
              -t    Wall time (HH:MM:SS).
              -m    Memory (e.g. 10G).
              -n    Number of CPUs.
              -N    Number of nodes.


DEFAULTS
========

  -subbasis:     No substitutions are done by default.
  -subpot:       No substitutions are done by default.
  -modpot:       No potentials are modified by default.
  -basisdb:      $OLCAO_DATA/atomicBDB
  -potdb:        $OLCAO_DATA/atomicPDB
  -scfkp:        One general kpoint is used for all SCF calculations.
  -pscfkp:       One general kpoint is used for all post-SCF calculations.
  -kp:           One general kpoint is used for all SCF and post-SCF calcs.
  -scfkpd:       No density-based SCF k-points by default.
  -pscfkpd:      No density-based post-SCF k-points by default.
  -kpd:          No density-based k-points by default.
  -scfkpint:     Gaussian/histogram (0) by default.
  -pscfkpint:    Gaussian/histogram (0) by default.
  -kpint:        Gaussian/histogram (0) by default.
  -kpshift:      Use the default shift for each lattice type.
  -printbz:      Use the first Brillouin zone and no scaling.
  -reduce:       Off by default.  Sub-parameter defaults:
                   level=2, thick=0.1 A, cutoff=4.0 A, operand=species,
                   tolerance=0.05, selection=0 (all atoms).
  -target:       Off by default.  Sub-parameter defaults:
                   zone=in, operand=species, relation=diff, radius=3.50 A.
                   The location is undefined.
  -block:        Off by default.  Sub-parameter defaults:
                   zone=in, operand=species, relation=diff.
                   The area is undefined.
  -xanes:        Off by default (one set of input files).
                   Default radius=3.50 A.  Atoms chosen 1 from each species.
  -xccode:       100 (Wigner interpolation method).
  -xcmesh:       numvect=100, weights=0.5/0.5, in_samp=0.1, out_samp=3.5,
                   spacing=0.8.
  -rel:          Non-relativistic by default.
  -statefactor:  2.5
  -pdb:          No PDB file generated.
  -cif:          No CIF file generated.
  -basisVis:     No visualization files generated.
  -nocore:       Core orbitals are orthogonalized out (default behavior).
  -slurm:        partition=rulisp-lab, account=rulisp-lab,
                   time=00:60:00, memory=10G, cpus=1, nodes=1.

Defaults are given in ./makeinputrc.py or $OLCAO_RC/makeinputrc.py.
"""

        parser = ap.ArgumentParser(
            prog=prog_name,
            formatter_class=ap.RawDescriptionHelpFormatter,
            description=description_text,
            epilog=epilog_text,
        )

        self._add_parser_arguments(parser)
        return parser.parse_args()

    def _add_parser_arguments(self, parser):
        """Register every CLI flag with argparse.

        Each argument's help string provides a concise summary.  For the full
        explanation of every option (including sub-options, interactions, and
        physics context), see the epilog printed by ``-h`` or the module
        docstring at the top of this file.
        """

        # ---- Database paths ----
        # The -basisdb option allows the user to supply the location of
        # another atomic basis functions database.  The default value is
        # $OLCAO_DATA/atomicBDB.
        parser.add_argument("-basisdb", dest="basisdb", type=str,
                            default=None,
                            help="Path to alternate atomic basis function "
                                 "database.  Default: $OLCAO_DATA/atomicBDB.")
        # The -potdb option allows the user to supply the location of
        # another potential function database.  The default value is
        # $OLCAO_DATA/atomicPDB.
        parser.add_argument("-potdb", dest="potdb", type=str,
                            default=None,
                            help="Path to alternate atomic potential function "
                                 "database.  Default: $OLCAO_DATA/atomicPDB.")

        # ---- Potential modification ----
        # The -modpot option adjusts the potential for a single element,
        # re-defining it with a new min/max alpha range and term count.
        # Only one element can be modified per invocation.
        parser.add_argument("-modpot", dest="modpot", nargs=4,
                            metavar=("ELEM", "MIN", "MAX", "NTERMS"),
                            default=None,
                            help="Modify potential for ELEM with new Gaussian "
                                 "alpha range [MIN, MAX] and NTERMS terms.  "
                                 "Only one element can be modified.")

        # ---- Basis / potential substitution (repeatable) ----
        # The -subbasis option substitutes an alternate basis definition.
        # If OUT is an element name (e.g. si) all atoms of that element are
        # substituted; if OUT is a species tag (e.g. si1) only that species
        # is substituted.  Example: "-subbasis b 2" uses contract2.dat
        # instead of contract1.dat for all boron.  Repeatable.
        parser.add_argument("-subbasis", dest="subbasis", nargs=2,
                            metavar=("OUT", "IN"), action="append",
                            default=None,
                            help="Substitute basis: replace OUT with IN.  "
                                 "OUT can be an element (si) or species "
                                 "(si1).  Repeatable.")
        # The -subpot option behaves identically to -subbasis but
        # substitutes the potential instead of the basis.  Repeatable.
        parser.add_argument("-subpot", dest="subpot", nargs=2,
                            metavar=("OUT", "IN"), action="append",
                            default=None,
                            help="Substitute potential: replace OUT with IN.  "
                                 "Same semantics as -subbasis.  Repeatable.")

        # ---- K-points ----
        # The -scfkp option specifies the k-point mesh for SCF calculations.
        # A B C are the mesh dimensions along the a, b, c lattice directions.
        # If all three are 1, the gamma point is assumed (faster: real-only
        # integral matrices).  If not given, 1 general kpoint is used.
        parser.add_argument("-scfkp", dest="scfkp", nargs=3, type=int,
                            metavar=("A", "B", "C"), default=None,
                            help="SCF k-point mesh dimensions.  1 1 1 = "
                                 "gamma point (faster, real matrices).")
        # The -pscfkp option is the same as -scfkp but applies to post-SCF
        # calculations (band, bond, dos, optc).  SYBD uses its own path.
        parser.add_argument("-pscfkp", dest="pscfkp", nargs=3, type=int,
                            metavar=("A", "B", "C"), default=None,
                            help="Post-SCF k-point mesh dimensions.  Same "
                                 "semantics as -scfkp but for band/bond/"
                                 "dos/optc.")
        # The -kp option applies the same mesh to both SCF and post-SCF.
        parser.add_argument("-kp", dest="kp", nargs=3, type=int,
                            metavar=("A", "B", "C"), default=None,
                            help="K-point mesh for both SCF and post-SCF.")
        # The -scfkpd option specifies a k-point volume
        # density for the SCF calculation.  The value is
        # the number of kpoints per unit reciprocal-space
        # volume (Bohr^-3).  Total kpoints = D * V_BZ.
        # OLCAO distributes them as uniformly as possible
        # across the three axes at runtime.  This bypasses
        # the makeKPoints program entirely.
        parser.add_argument(
            "-scfkpd", dest="scfkpd", type=float,
            metavar="D", default=None,
            help="SCF k-point volume density "
                 "(kpoints/Bohr^-3).  Total kpoints "
                 "= D * V_BZ.  Bypasses makeKPoints.")
        # The -pscfkpd option is the same but for post-SCF.
        parser.add_argument(
            "-pscfkpd", dest="pscfkpd", type=float,
            metavar="D", default=None,
            help="Post-SCF k-point volume density.  "
                 "Same semantics as -scfkpd but for "
                 "band/bond/dos/optc.")
        # The -kpd option sets the same density for both
        # SCF and post-SCF.  Density options are mutually
        # exclusive with explicit mesh options.
        parser.add_argument(
            "-kpd", dest="kpd", type=float,
            metavar="D", default=None,
            help="K-point volume density for both "
                 "SCF and post-SCF (kpoints/Bohr^-3)."
                 "  Mutually exclusive with "
                 "-kp/-scfkp/-pscfkp.")
        # The -scfkpint option selects the k-point
        # integration method for the SCF calculation.
        # 0 = Gaussian/histogram (default), 1 = LAT
        # (linear analytic tetrahedron).  Higher integers
        # are reserved for future methods.  This value is
        # written into the KPOINT_INTG_CODE field of the
        # k-point input file.
        parser.add_argument(
            "-scfkpint", dest="scfkpint", type=int,
            metavar="M", default=None,
            help="SCF k-point integration method.  "
                 "0 = Gaussian (default), 1 = LAT.")
        # The -pscfkpint option is the same but for
        # post-SCF calculations.
        parser.add_argument(
            "-pscfkpint", dest="pscfkpint", type=int,
            metavar="M", default=None,
            help="Post-SCF k-point integration method."
                 "  Same semantics as -scfkpint but "
                 "for band/bond/dos/optc.")
        # The -kpint option sets the same integration
        # method for both SCF and post-SCF.
        parser.add_argument(
            "-kpint", dest="kpint", type=int,
            metavar="M", default=None,
            help="K-point integration method for both "
                 "SCF and post-SCF.  0 = Gaussian "
                 "(default), 1 = LAT.")
        # The -kpshift option forces a specific fractional shift on the
        # k-point mesh.
        parser.add_argument("-kpshift", dest="kpshift", nargs=3, type=float,
                            metavar=("A", "B", "C"), default=None,
                            help="Fractional shift applied to the k-point "
                                 "mesh in a, b, c.")
        # The -printbz option writes a Brillouin zone description file.
        # BZ is the zone index; SCALE is a visualization scale factor.
        parser.add_argument("-printbz", dest="printbz", nargs=2,
                            metavar=("BZ", "SCALE"), default=None,
                            help="Write a BZ description file.  BZ = zone "
                                 "index, SCALE = visualization scale factor.")

        # ---- Exchange correlation ----
        # The -xccode option sets the exchange-correlation functional.
        # See $OLCAO_DIR/share/xc_code.dat for the list of codes.
        # Default is 100 (Wigner interpolation).
        parser.add_argument("-xccode", dest="xccode", type=int, default=None,
                            help="Exchange-correlation code (see xc_code.dat)."
                                 "  Default: 100 (Wigner interpolation).")
        # The -xcmesh flag enables the XC mesh sub-options below.  The mesh
        # is spherical and atom-centered, used for sampling the charge
        # density for exchange-correlation evaluation.
        parser.add_argument("-xcmesh", dest="xcmesh", action="store_true",
                            default=False,
                            help="Enable XC real-space mesh construction.  "
                                 "The mesh is spherical and atom-centered.")
        # The -numvect sub-option sets the number of radial sample vector
        # directions in the XC mesh.  Default: 100.
        parser.add_argument("-numvect", dest="numvect", type=int, default=None,
                            help="Number of radial sample vectors for the XC "
                                 "mesh.  Default: 100.")
        # The -weight sub-option sets the weighting of XC mesh points
        # inside (IN) and outside (OUT) a radial cutoff.  Default: 0.5, 0.5.
        parser.add_argument("-weight", dest="weight", nargs=2, type=float,
                            metavar=("IN", "OUT"), default=None,
                            help="XC mesh weights inside/outside radial "
                                 "cutoff.  Default: 0.5 0.5.")
        # The -samp sub-option defines the distribution of XC mesh points
        # along the vector directions: in-sampling, out-sampling, spacing.
        # Defaults: 0.1, 3.5, 0.8.
        parser.add_argument("-samp", dest="samp", nargs=3, type=float,
                            metavar=("IN", "OUT", "SPACING"), default=None,
                            help="XC mesh sampling: in-rate, out-rate, "
                                 "spacing.  Default: 0.1 3.5 0.8.")

        # ---- Target / Block / Reduce (repeatable) ----
        # These three options implement the atom-grouping schemes described
        # in the module docstring.  Each can be given multiple times on the
        # command line; they are applied in the order given.  Sub-options
        # are collected as a token list and parsed in reconcile().

        # The -target option groups atoms based on proximity to a specific
        # point (atom site, xyz, or abc coordinates).  Sub-options:
        #   -atom N        Target atom number N.
        #   -atxyz X Y Z   Target at Cartesian coordinates.
        #   -atabc A B C   Target at fractional coordinates.
        #   -sphere R      Radius of influence (default 3.50 A).
        #   -zone in|out   Consider atoms inside or outside the sphere.
        #   -operand X     Group by species, type, or reduce.
        #   -relate X      Make atoms alike or diff within the zone.
        parser.add_argument("-target", dest="target", action="append",
                            nargs=ap.REMAINDER, default=None,
                            help="Target grouping: group atoms near a point.  "
                                 "Sub-opts: -atom/-atxyz/-atabc, -sphere, "
                                 "-zone, -operand, -relate.  Repeatable.")
        # The -block option groups atoms within a 3D slab defined by
        # from/to values along each lattice direction.  Sub-options:
        #   -abc fromA toA fromB toB fromC toC
        #   -zone in|out
        #   -operand species|type
        #   -relate alike|diff
        # Letters 'a', 'b', 'c' can replace numeric "to" values to mean
        # the full lattice extent in that direction.
        parser.add_argument("-block", dest="block", action="append",
                            nargs=ap.REMAINDER, default=None,
                            help="Block grouping: group atoms in a 3D slab.  "
                                 "Sub-opts: -abc, -zone, -operand, -relate.  "
                                 "Repeatable.")
        # The -reduce option groups atoms by local-environment similarity
        # using concentric spherical shells.  Sub-options:
        #   -level N        Number of neighbor shells (default 2).
        #   -thick T        Shell thickness in Angstroms (default 0.1).
        #   -cutoff C       Max radius in Angstroms (default 4.0).
        #   -operand X      species (default) or type.
        #   -tolerance D    Distance match threshold in Ang (default 0.05).
        #   -selection S    Atom subset (0=all, default).
        # This method only makes species out of elements (not types out of
        # species).  It has no zone/relation parameters.
        parser.add_argument("-reduce", dest="reduce", action="append",
                            nargs=ap.REMAINDER, default=None,
                            help="Reduce grouping: group by local-environment "
                                 "similarity.  Sub-opts: -level, -thick, "
                                 "-cutoff, -operand, -tolerance, -selection.  "
                                 "Repeatable.")

        # ---- XANES ----
        # The -xanes option causes the script to generate one set of input
        # files for each target atom (one from each species by default, or
        # explicitly listed atoms).  Each set differs in that the target
        # atom has core orbitals included and nearby atoms get unique types.
        # Sub-options:
        #   -sphere R       Radius for unique-type assignment (default 3.50 A).
        #   -atom N [M..]   Specific target atom numbers (supports N..M ranges).
        parser.add_argument("-xanes", dest="xanes", nargs=ap.REMAINDER,
                            default=None,
                            help="Generate per-atom XANES input file sets.  "
                                 "Sub-opts: -sphere R, -atom N [M..].  "
                                 "Default radius 3.50 A.")

        # ---- SYBD path ----
        # The -sybdpath option specifies which reciprocal-space path to use
        # for the symmetric band structure (SYBD) calculation.  Valid names
        # are directories in $OLCAO_DATA/sybdDB.  Should be consistent with
        # the actual cell type.
        parser.add_argument("-sybdpath", dest="sybdpath", type=str,
                            default=None,
                            help="Name of the symmetric band path (from "
                                 "$OLCAO_DATA/sybdDB).  Must match the "
                                 "actual cell type.")

        # ---- Misc flags ----
        # The -rel option prepares input for a relativistic calculation.
        # (Not yet operational in the Fortran programs.)  Also increases
        # the default SCF iteration limit from 50 to 150.
        parser.add_argument("-rel", dest="rel", action="store_true",
                            default=False,
                            help="Prepare for relativistic calculation.  "
                                 "Also sets num_iter_main to 150.")
        # The -statefactor option controls how many states to calculate:
        # num_states = statefactor * num_valence_electrons.  Default: 2.5.
        parser.add_argument("-statefactor", dest="statefactor", type=float,
                            default=None,
                            help="State factor: num_states = factor * "
                                 "num_valence_electrons.  Default: 2.5.")
        # The -bfcut option sets the exponent cutoff for negligible
        # Gaussian basis function interactions.  Default: 16.
        parser.add_argument("-bfcut", dest="bfcut", type=float, default=None,
                            help="Basis function interaction cutoff exponent.  "
                                 "Default: 16.")
        # The -escort option sets the exponent cutoff for negligible
        # electrostatic interactions.  Default: 16.
        parser.add_argument("-escort", dest="escort", type=float, default=None,
                            help="Electrostatic interaction cutoff exponent.  "
                                 "Default: 16.")
        # The -pdb option generates a PDB (Protein Data Bank) crystal
        # structure file for visualization in external programs.
        parser.add_argument("-pdb", dest="pdb", action="store_true",
                            default=False,
                            help="Generate a PDB crystal structure file.")
        # The -cif option generates a CIF (Crystallographic Interchange
        # Format) structure file.  Preferred over PDB because it supports
        # greater numbers of species per atom.
        parser.add_argument("-cif", dest="cif_flag", action="store_true",
                            default=False,
                            help="Generate a CIF structure file (preferred "
                                 "over PDB; supports more species).")
        # The -basisVis option asks the contract program to produce files
        # in .inputTemp for visualizing the atomic orbital basis: numerical
        # radial data and POVRay scene files for each orbital type.
        parser.add_argument("-basisVis", dest="basis_vis", action="store_true",
                            default=False,
                            help="Generate basis visualization files: "
                                 "numerical radial data + POVRay scenes "
                                 "in .inputTemp/.")
        # The -emu option generates EMU configuration files and skips
        # the standard OLCAO output.
        parser.add_argument("-emu", dest="emu", action="store_true",
                            default=False,
                            help="Generate EMU configuration files.")
        # The -nocore option includes all core orbitals in the valence
        # (no orthogonalization).  This allows computation of core energy
        # eigenvalues.  Without this flag, "core" basis functions are
        # orthogonalized out of the eigenvalue problem.
        parser.add_argument("-nocore", dest="nocore", action="store_true",
                            default=False,
                            help="Include all core orbitals in the valence "
                                 "(no orthogonalization).  Allows computing "
                                 "core eigenvalues.")

        # ---- Slurm ----
        # The -slurm option specifies parameters for the slurm job
        # submission script.  Sub-options:
        #   -p PARTITION    Partition name (default: rulisp-lab).
        #   -a ACCOUNT      Account to charge (default: rulisp-lab).
        #   -t TIME         Wall time as HH:MM:SS (default: 00:60:00).
        #   -m MEMORY       Memory request (default: 10G).
        #   -n CPUS         Number of CPUs (default: 1).
        #   -N NODES        Number of nodes (default: 1).
        parser.add_argument("-slurm", dest="slurm", nargs=ap.REMAINDER,
                            default=None,
                            help="Slurm submission parameters.  Sub-opts: "
                                 "-p partition, -a account, -t time, "
                                 "-m memory, -n cpus, -N nodes.")

    # ------------------------------------------------------------------
    # Reconcile CLI with rc defaults
    # ------------------------------------------------------------------

    def reconcile(self, args):
        """Merge parsed CLI arguments into the settings object."""

        # Database paths.
        if args.basisdb is not None:
            self.atomic_bdb = args.basisdb
        if args.potdb is not None:
            self.atomic_pdb = args.potdb

        # Potential modification.
        if args.modpot is not None:
            self.mod_pot = 1
            self.mod_element_name = args.modpot[0].lower()
            self.min_mod_term = float(args.modpot[1])
            self.max_mod_term = float(args.modpot[2])
            self.num_mod_terms = int(args.modpot[3])

        # Basis / potential substitutions.
        if args.subbasis is not None:
            for out_val, in_val in args.subbasis:
                self.basis_sub_out.append(out_val)
                self.basis_sub_in.append(in_val)
        if args.subpot is not None:
            for out_val, in_val in args.subpot:
                self.pot_sub_out.append(out_val)
                self.pot_sub_in.append(in_val)

        # K-points -- density mode vs. explicit mesh mode.
        # Check whether any density option was given.
        has_density = (args.kpd is not None
                       or args.scfkpd is not None
                       or args.pscfkpd is not None)
        has_mesh = (args.kp is not None
                    or args.scfkp is not None
                    or args.pscfkp is not None)

        if has_density:
            # Density mode is all-or-nothing.  If any
            # density flag is set, both groups use density.
            self.use_kp_density = True
            if has_mesh:
                print("WARNING: Both density (-kpd/"
                      "-scfkpd/-pscfkpd) and explicit "
                      "mesh (-kp/-scfkp/-pscfkp) options "
                      "were given.  Density mode wins; "
                      "explicit mesh options are ignored.")

            # Default density is 1 for any unset group.
            self.kp_density[1] = 1.0
            self.kp_density[2] = 1.0

            # Apply the user's density values.
            if args.kpd is not None:
                self.kp_density[1] = args.kpd
                self.kp_density[2] = args.kpd
            if args.scfkpd is not None:
                self.kp_density[1] = args.scfkpd
            if args.pscfkpd is not None:
                self.kp_density[2] = args.pscfkpd

            # Update kp_note to reflect density mode.
            self.kp_note[1] = "(Density)"
            self.kp_note[2] = "(Density)"

            # Warn if -printbz was requested; it requires
            # the makeKPoints executable which is bypassed
            # in density mode.
            if args.printbz is not None:
                print("NOTE: -printbz is being skipped. "
                      "Brillouin zone visualization only "
                      "works with an explicit k-point "
                      "mesh (-kp/-scfkp/-pscfkp).")
        else:
            # Explicit mesh mode (original behavior).
            if args.kp is not None:
                self.kp_mesh_scf = list(args.kp)
                self.kp_mesh_pscf = list(args.kp)
                if args.kp == [1, 1, 1]:
                    self.kp_note[1] = "(Gamma)"
                    self.kp_note[2] = "(Gamma)"
            if args.scfkp is not None:
                self.kp_mesh_scf = list(args.scfkp)
                if args.scfkp == [1, 1, 1]:
                    self.kp_note[1] = "(Gamma)"
            if args.pscfkp is not None:
                self.kp_mesh_pscf = list(args.pscfkp)
                if args.pscfkp == [1, 1, 1]:
                    self.kp_note[2] = "(Gamma)"
            if args.printbz is not None:
                self.print_bz = int(args.printbz[0])
                self.scale_factor = float(
                    args.printbz[1])

        # K-point integration method.  Apply the combined
        # flag first, then let per-group flags override.
        if args.kpint is not None:
            self.kp_intg_code[1] = args.kpint
            self.kp_intg_code[2] = args.kpint
        if args.scfkpint is not None:
            self.kp_intg_code[1] = args.scfkpint
        if args.pscfkpint is not None:
            self.kp_intg_code[2] = args.pscfkpint

        # Non-Gaussian integration methods (e.g. LAT)
        # require density mode because the makeKPoints
        # Fortran program hardcodes KPOINT_INTG_CODE = 0.
        # Once mesh generation moves into OLCAO proper
        # this restriction will be lifted.
        has_nonzero_intg = (
            self.kp_intg_code[1] != 0
            or self.kp_intg_code[2] != 0)
        if has_nonzero_intg and not self.use_kp_density:
            intg_names = {0: "Gaussian", 1: "LAT"}
            scf_name = intg_names.get(
                self.kp_intg_code[1],
                str(self.kp_intg_code[1]))
            pscf_name = intg_names.get(
                self.kp_intg_code[2],
                str(self.kp_intg_code[2]))
            print(
                f"ERROR: -kpint/-scfkpint/-pscfkpint "
                f"with a non-Gaussian method "
                f"(SCF={scf_name}, PSCF={pscf_name}) "
                f"requires density mode (-kpd/-scfkpd/"
                f"-pscfkpd).  The legacy makeKPoints "
                f"program does not support alternative "
                f"integration methods.  Use density "
                f"mode or omit the -kpint flag.")
            sys.exit(1)

        # Shift applies to both mesh and density modes.
        if args.kpshift is not None:
            self.kp_shift = (
                f"{args.kpshift[0]} "
                f"{args.kpshift[1]} "
                f"{args.kpshift[2]}")

        # Exchange correlation.
        if args.xccode is not None:
            self.xc_code = args.xccode
        if args.numvect is not None:
            self.num_samp_vectors = args.numvect
        if args.weight is not None:
            self.xc_in_weight = args.weight[0]
            self.xc_out_weight = args.weight[1]
        if args.samp is not None:
            self.xc_in_samp = args.samp[0]
            self.xc_out_samp = args.samp[1]
            self.xc_spacing_samp = args.samp[2]

        # Grouping methods: reduce.
        if args.reduce is not None:
            for sub_args in args.reduce:
                self._parse_reduce(sub_args)

        # Grouping methods: target.
        if args.target is not None:
            for sub_args in args.target:
                self._parse_target(sub_args)

        # Grouping methods: block.
        if args.block is not None:
            for sub_args in args.block:
                self._parse_block(sub_args)

        # XANES.
        if args.xanes is not None:
            self.xanes = 1
            self._parse_xanes(args.xanes)

        # SYBD path.
        if args.sybdpath is not None:
            self.sybd_path = args.sybdpath.lower()

        # Misc flags.
        if args.rel:
            self.rel = 1
            self.num_iter_main = 150
        if args.statefactor is not None:
            self.state_factor = args.statefactor
        if args.bfcut is not None:
            self.bf_cutoff = args.bfcut
        if args.escort is not None:
            self.es_cutoff = args.escort
        if args.pdb:
            self.pdb = 1
        if args.cif_flag:
            self.cif = 1
        if args.basis_vis:
            self.do_basis_vis = 1
        if args.emu:
            self.emu = 1
        if args.nocore:
            self.no_core = 1

        # Slurm.
        if args.slurm is not None:
            self._parse_slurm(args.slurm)

    # ------------------------------------------------------------------
    # Sub-option parsers for compound flags
    # ------------------------------------------------------------------

    def _parse_reduce(self, tokens):
        """Parse sub-options for one -reduce invocation."""
        from makeinputrc import parameters_and_defaults
        rc = parameters_and_defaults()
        r = {
            "level":     rc["reduce_level"],
            "thick":     rc["reduce_thick"],
            "cutoff":    rc["reduce_cutoff"],
            "op":        rc["reduce_op"],
            "tolerance": rc["reduce_tolerance"],
            "selection": rc["reduce_selection"],
        }
        i = 0
        while i < len(tokens):
            tag = tokens[i]
            if tag == "-level":
                i += 1; r["level"] = int(tokens[i])
            elif tag == "-thick":
                i += 1; r["thick"] = float(tokens[i])
            elif tag == "-cutoff":
                i += 1; r["cutoff"] = float(tokens[i])
            elif tag == "-operand":
                i += 1; r["op"] = tokens[i].lower()
            elif tag == "-tolerance":
                i += 1; r["tolerance"] = float(tokens[i])
            elif tag == "-selection":
                i += 1; r["selection"] = int(tokens[i])
            else:
                break
            i += 1
        self.reduces.append(r)
        self.methods.append(("reduce", len(self.reduces) - 1))

    def _parse_target(self, tokens):
        """Parse sub-options for one -target invocation."""
        from makeinputrc import parameters_and_defaults
        rc = parameters_and_defaults()
        t = {
            "radius":   rc["target_radius"],
            "zone":     rc["target_zone"],
            "op":       rc["target_op"],
            "relation": rc["target_relation"],
            "loc":      None,
            "loc_type": None,   # 1=atom, 2=xyz, 3=abc
        }
        i = 0
        while i < len(tokens):
            tag = tokens[i]
            if tag == "-atom":
                i += 1; t["loc"] = int(tokens[i]); t["loc_type"] = 1
            elif tag == "-atxyz":
                t["loc"] = [float(tokens[i+1]), float(tokens[i+2]),
                            float(tokens[i+3])]
                t["loc_type"] = 2; i += 3
            elif tag == "-atabc":
                t["loc"] = [float(tokens[i+1]), float(tokens[i+2]),
                            float(tokens[i+3])]
                t["loc_type"] = 3; i += 3
            elif tag == "-sphere":
                i += 1; t["radius"] = float(tokens[i])
            elif tag == "-zone":
                i += 1; t["zone"] = tokens[i].lower()
            elif tag == "-operand":
                i += 1; t["op"] = tokens[i].lower()
            elif tag == "-relate":
                i += 1; t["relation"] = tokens[i].lower()
            else:
                break
            i += 1
        self.targets.append(t)
        self.methods.append(("target", len(self.targets) - 1))

    def _parse_block(self, tokens):
        """Parse sub-options for one -block invocation."""
        from makeinputrc import parameters_and_defaults
        rc = parameters_and_defaults()
        b = {
            "zone":     rc["block_zone"],
            "op":       rc["block_op"],
            "relation": rc["block_relation"],
            "borders":  None,   # [[from_a,to_a],[from_b,to_b],[from_c,to_c]]
        }
        i = 0
        while i < len(tokens):
            tag = tokens[i]
            if tag == "-abc":
                b["borders"] = [
                    [tokens[i+1], tokens[i+2]],
                    [tokens[i+3], tokens[i+4]],
                    [tokens[i+5], tokens[i+6]],
                ]
                i += 6
            elif tag == "-zone":
                i += 1; b["zone"] = tokens[i].lower()
            elif tag == "-operand":
                i += 1; b["op"] = tokens[i].lower()
            elif tag == "-relate":
                i += 1; b["relation"] = tokens[i].lower()
            else:
                break
            i += 1
        self.blocks.append(b)
        self.methods.append(("block", len(self.blocks) - 1))

    def _parse_xanes(self, tokens):
        """Parse sub-options for -xanes."""
        i = 0
        while i < len(tokens):
            tag = tokens[i]
            if tag == "-sphere":
                i += 1; self.xanes_radius = float(tokens[i])
            elif tag == "-atom":
                i += 1
                while i < len(tokens) and not tokens[i].startswith("-"):
                    tok = tokens[i]
                    if ".." in tok:
                        lo, hi = tok.split("..")
                        self.xanes_atoms.extend(range(int(lo), int(hi) + 1))
                    else:
                        self.xanes_atoms.append(int(tok))
                    i += 1
                continue  # skip the i += 1 at the bottom
            else:
                break
            i += 1

    def _parse_slurm(self, tokens):
        """Parse sub-options for -slurm."""
        i = 0
        while i < len(tokens):
            tag = tokens[i]
            if tag == "-p":
                i += 1; self.partition = tokens[i]
            elif tag == "-a":
                i += 1; self.account = tokens[i]
            elif tag == "-t":
                i += 1; self.time = tokens[i]
            elif tag == "-m":
                i += 1; self.memory = tokens[i]
            elif tag == "-n":
                i += 1; self.cpus = int(tokens[i])
            elif tag == "-N":
                i += 1; self.nodes = int(tokens[i])
            else:
                break
            i += 1

    # ------------------------------------------------------------------
    # Record keeping
    # ------------------------------------------------------------------

    def record_clp(self):
        """Append the command line used to the ``command`` file."""
        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime("%b. %d, %Y: %H:%M:%S")
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write("Cmnd:")
            for arg in sys.argv:
                cmd.write(f" {arg}")
            cmd.write("\n\n")


# ---------------------------------------------------------------------------
# File-name constants
# ---------------------------------------------------------------------------

# Directory names.
INPUT_TEMP = ".inputTemp"
INPUTS_DIR = "inputs"

# File names for OLCAO input generation.
OLCAO_SKL      = "olcao.skl"
OLCAO_FMI      = "olcao.fract-mi"
OLCAO_CMI      = "olcao.cart-mi"
OLCAO_DAT      = "olcao.dat"
POTENTIAL      = "scfV.dat"
KP_SCF         = "kp-scf.dat"
KP_PSCF        = "kp-pscf.dat"
KP_IN_FILE     = "kpSpecs.dat"
KP_OUT_FILE    = "kpSpecs.out"
STRUCTURE      = "structure.dat"
EMU_CONF       = "emu.conf"
PSO_CONF       = "pso.conf"

# Tracking / submission files.
PBS_SUB        = "pbs"
LSF_SUB        = "lsf"
SLURM_SUB      = "slurm"
BASH_SUB       = "bash"
COMMAND_FILE   = "command"
SUMMARY_FILE   = "summary"
MEMORY_FILE    = "memory"
DAT_SKL_MAP    = "datSkl.map"
BATCH_LIST     = "batchList"

# File extensions for structure visualisation.
PDB_EXT        = ".pdb"
CIF_EXT        = ".cif"

# Size of one OLCAO double-precision number in bytes (set by OLCAO_DOUBLE
# environment variable, default 16 = Fortran COMPLEX*16).
OLCAO_DOUBLE   = int(os.getenv("OLCAO_DOUBLE", "16"))

# XANES edge tags (indexed 0..16).
EDGE_TAGS = [
    "",    "1s_", "2s_", "2p_", "3s_", "3p_", "4s_", "3d_",
    "4p_", "5s_", "4d_", "5p_", "6s_", "4f_", "5d_", "6p_", "7s_",
]

# XANES initial states table: init_xanes_states[n][l] = (first, last).
INIT_XANES_STATES = {
    1: {0: (1, 1)},                                          # 1s
    2: {0: (2, 2),   1: (3, 5)},                             # 2s, 2p
    3: {0: (6, 6),   1: (7, 9),   2: (10, 14)},             # 3s, 3p, 3d
    4: {0: (15, 15), 1: (16, 18), 2: (19, 23), 3: (24, 30)},# 4s-4f
    5: {0: (31, 31), 1: (32, 34), 2: (35, 39)},             # 5s, 5p, 5d
    6: {0: (40, 40), 1: (41, 43)},                           # 6s, 6p
    7: {0: (44, 44)},                                        # 7s
}


# ---------------------------------------------------------------------------
# Environment setup
# ---------------------------------------------------------------------------

def remove_old_input():
    """Delete stale input files from a previous run."""
    for f in [OLCAO_DAT, STRUCTURE, KP_SCF, KP_PSCF]:
        if os.path.exists(f):
            os.remove(f)
    for tag in EDGE_TAGS:
        fname = f"{tag}{POTENTIAL}"
        if os.path.exists(fname):
            os.remove(fname)


def setup_environment(settings):
    """Clean previous outputs, set up directories, resolve database paths."""
    print("      Cleaning up Previous Inputs")
    remove_old_input()

    print("      Initializing Directories")
    os.makedirs(INPUT_TEMP, exist_ok=True)

    # Resolve database paths (may already be set by CLI).
    olcao_data = os.getenv("OLCAO_DATA", "")
    if not hasattr(settings, "atomic_bdb"):
        settings.atomic_bdb = os.path.join(olcao_data, "atomicBDB")
    if not hasattr(settings, "atomic_pdb"):
        settings.atomic_pdb = os.path.join(olcao_data, "atomicPDB")
    settings.sybd_db = os.path.join(olcao_data, "sybdDB")
    settings.space_db = os.path.join(olcao_data, "spaceDB")


# ---------------------------------------------------------------------------
# Main workflow steps -- stubs to be filled in progressively
# ---------------------------------------------------------------------------

def initialize_cell(settings, sc):
    """Read olcao.skl, apply space group and supercell, init species.

    This subroutine reads the olcao.skl skeleton input file and initialises
    all the parameters that define the structure.  The StructureControl
    method ``read_input_file`` handles the heavy lifting: parsing the
    skeleton, applying the space group, applying the supercell, creating
    the element list, creating the species data, mapping element numbers,
    and computing implicit information (potential alphas, valence states).

    After the structure is loaded we:
      1. Copy key structural data into the settings object for convenient
         access later (element/species/type IDs, lattice data, etc.).
      2. Initialise the type-ID arrays so that every atom starts with
         type 1 within its species (types may be reassigned later by the
         grouping methods: target, block, reduce, or XANES).
      3. Determine the default symmetric band structure (SYBD) path from
         the space group and lattice type, following the conventions of
         Setyawan and Curtarolo, Comp. Mat. Sci. 49, pp. 299-312 (2010).
      4. Set the default k-point mesh shift based on the lattice type.
      5. Move the space-group scratch files (sginput, sgoutput) into the
         .inputTemp directory so they are preserved with the other
         intermediate files.

    Parameters
    ----------
    settings : ScriptSettings
        The command-line / rc-file settings object.
    sc : StructureControl
        The structure control instance (will be populated by this call).
    """
    from structure_control import PI, BOHR_RAD

    # A small tolerance for floating-point lattice comparisons.
    epsilon = 1e-5

    # Read olcao.skl and request to use the species defined in that file.
    # The use_file_species=True flag tells create_species_data to honour
    # the species tags in the skeleton (e.g. si1, si2) rather than
    # collapsing all atoms of each element into a single species.  This
    # preserves the user's explicit type assignments, which may then be
    # further refined by the grouping methods (target, block, reduce).
    sc.read_input_file(OLCAO_SKL, use_file_species=True)

    # ------------------------------------------------------------------
    # Transfer structural data into settings for convenient access.
    # In the Perl version these were global variables; here they live on
    # the settings object so that every downstream function can reach
    # them via the same ``settings`` reference that carries the CLI params.
    # ------------------------------------------------------------------

    # Scalar lattice / space-group metadata.
    settings.num_atoms         = sc.num_atoms
    settings.num_elements      = sc.num_elements
    settings.space_group_name  = sc.space_group_name
    settings.space_group_num   = sc.space_group_num
    settings.space_group_sub_num = sc.space_group_sub_num
    settings.lattice_type      = sc.lattice_type
    settings.do_full_cell      = sc.do_full_cell
    settings.min_pot_alpha     = sc.min_pot_alpha
    settings.max_num_pot_alphas = sc.max_num_pot_alphas
    settings.max_num_atom_alphas = sc.max_num_atom_alphas
    settings.max_num_vale_states = sc.max_num_vale_states

    # Per-atom arrays (1-indexed, index 0 = None).
    # We make copies so that the grouping methods can modify the local
    # versions (especially type IDs) without touching StructureControl.
    settings.atom_element_name = list(sc.atom_element_name)
    settings.atom_element_id   = list(sc.atom_element_id)
    settings.atom_species_id   = list(sc.atom_species_id)

    # Element list (1-indexed).
    settings.element_list = list(sc.element_list)

    # Species counts per element (1-indexed).
    settings.num_species = list(sc.num_species)

    # Initialise type IDs.  Every atom starts as type 1 within its
    # species.  atom_type_id is indexed [file_set][atom] where
    # file_set 0 is the base (non-XANES) configuration.  For XANES,
    # additional file sets will be appended later.
    settings.atom_type_id = [[None]]  # [0][0] placeholder
    for atom in range(1, settings.num_atoms + 1):
        settings.atom_type_id[0].append(1)

    # Initialise type counts per species per element.
    # num_types[file_set][element][species] — all start at 1.
    settings.num_types = [[None]]  # [0][0] placeholder
    for element in range(1, settings.num_elements + 1):
        elem_types = [None]
        for species in range(1, settings.num_species[element] + 1):
            elem_types.append(1)
        settings.num_types[0].append(elem_types)

    # ------------------------------------------------------------------
    # Determine the symmetric band structure (SYBD) path.
    # ------------------------------------------------------------------
    # The high-symmetry points of the first Brillouin zone (BZ) are given
    # in the sybd database in fractional coordinates of the reciprocal
    # space lattice vectors of the appropriate lattice type (conventional
    # or primitive).
    #
    # If an fcc cubic system is given and the primitive cell is used, then
    # the reciprocal lattice of the primitive cell is also used and so the
    # high symmetry points must come from the fcc reciprocal primitive
    # lattice.  If an fcc cubic system is given and the conventional
    # (full) cell is used, then the reciprocal lattice of the conventional
    # (full) cell is also used and so the high symmetry points must come
    # from the simple cubic reciprocal primitive lattice.  This holds for
    # the other crystal systems too.
    #
    # There are two cases: (1) no command line -sybdpath was given, so we
    # auto-detect from the space group and lattice; (2) the user gave an
    # explicit -sybdpath, so we use that directly.

    if settings.sybd_path == "":
        _auto_sybd_path(settings, sc, epsilon)

    # ------------------------------------------------------------------
    # Default k-point mesh shift (if not set on the command line).
    # ------------------------------------------------------------------
    if settings.kp_shift == "-1 -1 -1":
        _auto_kp_shift(settings)

    # ------------------------------------------------------------------
    # Move the space-group scratch files into .inputTemp.
    # ------------------------------------------------------------------
    _safe_move("sginput", os.path.join(INPUT_TEMP, "sginput"))
    _safe_move("sgoutput", os.path.join(INPUT_TEMP, "sgoutput"))


def _safe_move(src, dst):
    """Move ``src`` to ``dst`` if ``src`` exists; otherwise do nothing.

    If ``src`` does not exist but ``dst`` already does, print a notice
    and leave ``dst`` untouched (matches Perl ``move`` behaviour).
    """
    if not os.path.exists(src):
        if os.path.exists(dst):
            print(f"      {src} not found but {dst} exists.  No move done.")
        return
    shutil.move(src, dst)


def _auto_sybd_path(settings, sc, epsilon):
    """Determine the SYBD path automatically from the space group.

    This implements the path-selection logic described by Setyawan and
    Curtarolo, Comp. Mat. Sci. 49, pp. 299-312 (2010).  The selection
    depends on the space group number, the lattice type (P, F, I, C, R,
    H, A, B), whether the full cell is being used, and in some cases
    the reciprocal-lattice angles and lattice-vector magnitudes.

    Parameters
    ----------
    settings : ScriptSettings
        Settings object — ``sybd_path`` will be set on return.
    sc : StructureControl
        Populated structure control instance.
    epsilon : float
        Tolerance for floating-point comparisons (~1e-5).
    """
    pi = math.pi
    sg   = settings.space_group_num
    lt   = settings.lattice_type
    full = settings.do_full_cell

    # Reciprocal-lattice angles (1-indexed: [1]=alpha*, [2]=beta*, [3]=gamma*).
    ar = sc.angle_recip

    # Full-cell magnitudes and angles (1-indexed).
    fcm = sc.full_cell_mag
    fca = sc.full_cell_angle

    # Primitive-cell angles (1-indexed).
    ang = sc.angle

    # ------------------------------------------------------------------
    # Triclinic  (SG 1–2)
    # ------------------------------------------------------------------
    if sg <= 2:
        if (ar[1] > pi/2 and ar[2] > pi/2 and ar[3] > pi/2
                and ar[3] == min(ar[1], ar[2], ar[3])):
            settings.sybd_path = "tric1a"
        elif (ar[1] < pi/2 and ar[2] < pi/2 and ar[3] < pi/2
                and ar[3] == max(ar[1], ar[2], ar[3])):
            settings.sybd_path = "tric1b"
        elif (ar[1] > pi/2 and ar[2] > pi/2 and abs(ar[3] - pi/2) < epsilon):
            settings.sybd_path = "tric2a"
        elif (ar[1] < pi/2 and ar[2] < pi/2 and abs(ar[3] - pi/2) < epsilon):
            settings.sybd_path = "tric2b"
        else:
            print("      No SYBD path identifiable.")
            print("      Expecting form as in Setyawan and Curtarolo")
            print("      Comp. Mat. Sci., v49, pp299-312, (2010)")
            print("      Defaulting to use 'tric1a' as the path.")
            settings.sybd_path = "tric1a"

    # ------------------------------------------------------------------
    # Monoclinic  (SG 3–15)
    # ------------------------------------------------------------------
    elif sg <= 15:
        # Assign a sybd path for P, A, B, F, and I all as a "standard"
        # monoclinic path.  Only for the C-centered case will we make
        # some distinctions.  In the future this will need to be fixed
        # for the other settings.
        if lt in ("P", "A", "B", "I", "F") or full:
            settings.sybd_path = "mono"
        elif ar[3] > pi / 2:
            settings.sybd_path = "monoc1"
        elif abs(ar[3] - pi / 2) < epsilon:
            settings.sybd_path = "monoc2"
        else:
            # ar[3] < pi/2 — need to distinguish monoc3, monoc4, monoc5.
            factor = (fcm[2] * math.cos(fca[1]) / fcm[3]
                      + fcm[2]**2 * math.sin(fca[1])**2 / fcm[1]**2)
            # Out of order for easier programming and computation.
            if abs(factor - 1) < epsilon:
                settings.sybd_path = "monoc4"
            elif factor < 1:
                settings.sybd_path = "monoc3"
            else:
                settings.sybd_path = "monoc5"

    # ------------------------------------------------------------------
    # Orthorhombic  (SG 16–74)
    # ------------------------------------------------------------------
    elif sg <= 74:
        # Assign a sybd path for P, A, and B.
        if lt in ("P", "A", "B") or full:
            settings.sybd_path = "ortho"
        elif lt == "I":
            settings.sybd_path = "orthoi"
        elif lt == "C":
            settings.sybd_path = "orthoc"
        else:
            # Face-centered orthorhombic.
            factor1 = 1.0 / fcm[1]**2
            factor2 = 1.0 / fcm[2]**2 + 1.0 / fcm[3]**2
            # Out of order for easier programming and computation.
            if abs(factor1 - factor2) < epsilon:
                settings.sybd_path = "orthof3"
            elif factor1 > factor2:
                settings.sybd_path = "orthof1"
            else:
                settings.sybd_path = "orthof2"

    # ------------------------------------------------------------------
    # Tetragonal  (SG 75–142)
    # ------------------------------------------------------------------
    elif sg <= 142:
        # Assign sybd path for P, F, and C.
        if lt in ("P", "F", "C") or full:
            settings.sybd_path = "tet"
        elif lt == "I":
            # Body-centered tetragonal.
            if fcm[3] < fcm[1]:
                settings.sybd_path = "teti1"
            else:
                settings.sybd_path = "teti2"

    # ------------------------------------------------------------------
    # Trigonal  (SG 143–167)
    # ------------------------------------------------------------------
    elif sg <= 167:
        if lt in ("P", "H"):
            settings.sybd_path = "hex"
        elif lt == "R" and settings.space_group_sub_num == 2 and full:
            # Hexagonal cell made of 3 rhombohedral cells.
            settings.sybd_path = "hex"
        else:
            # Rhombohedral cell.  The selection between rhomb1 and rhomb2
            # depends on the primitive-cell angle alpha.
            #
            # In the Setyawan and Curtarolo 2010 paper they propose two
            # variations for the rhombohedral symmetric band structure
            # path.  The selection of which path to use is said to depend
            # on the value of the cell angle alpha.  In the paper, the
            # angle alpha is defined as the interaxial angle between the
            # b and c axes of the conventional cell.  For trigonal systems
            # the conventional cell is nominally hexagonal and so the
            # angle alpha is defined to always be 90 degrees.  But, in
            # the paper the alpha dependence states that one path should
            # be taken when alpha > 90 and another path taken when the
            # angle alpha < 90 degrees.  There is no "equal to 90 degrees"
            # option specified.  The combination of those facts indicates
            # that the angle alpha for the rhombohedral cell should be the
            # angle alpha of the primitive cell which follows alpha = beta
            # = gamma for the rhombohedral primitive cell.  (Also, for
            # that cell a = b = c.)  Therefore, even though the paper says
            # that alpha should be understood as the angle from the
            # conventional cell we should take it in this case as the
            # angle from the primitive cell.  Specifically, that is
            # angle[1].
            if ang[1] < pi / 2:
                settings.sybd_path = "rhomb1"  # includes prim cells of R*_b groups
            else:
                settings.sybd_path = "rhomb2"

    # ------------------------------------------------------------------
    # Hexagonal  (SG 168–194)
    # ------------------------------------------------------------------
    elif sg <= 194:
        settings.sybd_path = "hex"

    # ------------------------------------------------------------------
    # Cubic  (SG 195–230)
    # ------------------------------------------------------------------
    elif sg <= 230:
        if full or lt == "P":
            settings.sybd_path = "sc"
        elif lt == "F":
            settings.sybd_path = "fcc"
        elif lt == "I":
            settings.sybd_path = "bcc"
        else:
            settings.sybd_path = "sc"


def _auto_kp_shift(settings):
    """Set the default k-point mesh shift based on the SYBD path / lattice.

    The shift is only set here if the user did not provide an explicit
    -kpshift on the command line (indicated by the sentinel value
    "-1 -1 -1").  The shift values are lattice-type-dependent defaults.

    Parameters
    ----------
    settings : ScriptSettings
        Settings object — ``kp_shift`` and ``cell_name`` will be set.
    """
    sp = settings.sybd_path

    if "tric" in sp:
        settings.cell_name = "triclinic"
        settings.kp_shift  = "0.5 0.5 0.5"
    elif "mono" in sp:
        settings.cell_name = "monoclinic"
        settings.kp_shift  = "0.5 0.5 0.5"
    elif "ortho" in sp:
        settings.cell_name = "orthorhombic"
        settings.kp_shift  = "0.5 0.5 0.5"
    elif "tet" in sp:
        settings.cell_name = "tetragonal"
        settings.kp_shift  = "0.5 0.5 0.5"
    elif "rhomb" in sp:
        settings.cell_name = "rhombohedral"
        settings.kp_shift  = "0.5 0.5 0.5"
    elif sp == "hex":
        settings.cell_name = "hexagonal"
        settings.kp_shift  = "0.3333333333333333 0.3333333333333333 0.25"
    else:
        # Nothing but cubic systems from here on.
        settings.cell_name = "cubic"
        if sp == "sc":
            settings.kp_shift = "0.25 0.25 0.25"
        elif sp == "fcc":
            settings.kp_shift = "0.25 0.25 0.25"
        elif sp == "bcc":
            settings.kp_shift = "0.5 0.5 0.5"


def assign_group(settings, sc, group_type):
    """Assign atoms to a group (species or types) by applying the requested
    grouping methods in the order they appeared on the command line.

    This subroutine walks through ``settings.methods`` — the ordered list of
    (method_name, index) tuples built during command-line parsing — and for
    each one checks whether its operand matches the current ``group_type``.
    If it does, the corresponding grouping function is called.

    The three grouping methods are:

    - **reduce**: Groups atoms by local-environment similarity using
      concentric spherical shells.  Only operates on species (not types).
    - **target**: Groups atoms based on proximity to a specific point
      (atom site, xyz, or abc coordinates).
    - **block**: Groups atoms based on whether they fall inside or outside
      a 3D slab defined by from/to values along each lattice direction.

    After all applicable methods have been applied, ``relax_group`` is
    called to renumber the species or types, filling in any gaps left by
    the assignment process.  (For example, if Si3 was reassigned and no
    atoms remain with that species, the numbering is compressed so there
    are no holes.)

    Before any grouping can happen, the interaction limit distance is set
    (to the maximum of all target radii, reduce cutoffs, and the XANES
    radius), target coordinates are prepared, and the minimum-distance
    matrix is computed.

    Parameters
    ----------
    settings : ScriptSettings
        The settings object carrying CLI parameters and mutable state
        (atom IDs, species counts, type counts, etc.).
    sc : StructureControl
        The populated structure control instance.
    group_type : str
        Either ``"species"`` or ``"types"``.  This subroutine is called
        twice from main — once for each operand.
    """
    # ------------------------------------------------------------------
    # Compute the interaction limit distance.  This is the maximum of
    # all target radii, all reduce cutoffs, and the XANES radius.  It
    # is passed to StructureControl so that the minimum-distance matrix
    # covers a large enough range.
    # ------------------------------------------------------------------
    limit_dist = 4.0
    for t in settings.targets:
        if t["radius"] > limit_dist:
            limit_dist = t["radius"]
    for r in settings.reduces:
        if r["cutoff"] > limit_dist:
            limit_dist = r["cutoff"]
    if settings.xanes == 1:
        if settings.xanes_radius > limit_dist:
            limit_dist = settings.xanes_radius
    sc.set_limit_dist(limit_dist)

    # ------------------------------------------------------------------
    # Prepare target coordinate data (fractional abc and direct xyz) if
    # any targets were requested and we haven't done so yet.
    # ------------------------------------------------------------------
    if len(settings.targets) > 0 and not getattr(settings, '_targets_prepared', False):
        prepare_targets(settings, sc)

    # ------------------------------------------------------------------
    # Build the minimum-distance matrix if any grouping method needs it
    # and we haven't built it yet.
    # ------------------------------------------------------------------
    needs_matrix = (len(settings.targets) > 0 or settings.xanes == 1
                    or len(settings.reduces) > 0)
    if needs_matrix and not getattr(settings, '_min_dist_made', False):
        make_min_dist_matrices(settings, sc)

    # ------------------------------------------------------------------
    # Walk the methods list, dispatching to the appropriate grouper when
    # the operand matches the current group_type.
    # ------------------------------------------------------------------
    current_reduce = 0
    current_target = 0
    current_block  = 0
    relaxation_needed = False

    for method_name, method_idx in settings.methods:
        if method_name == "reduce":
            r = settings.reduces[method_idx]
            if r["op"] == group_type:
                group_reduce(settings, sc, method_idx)
                relaxation_needed = True
            current_reduce += 1
        elif method_name == "target":
            t = settings.targets[method_idx]
            if t["op"] == group_type:
                group_target(settings, sc, method_idx)
                relaxation_needed = True
            current_target += 1
        elif method_name == "block":
            b = settings.blocks[method_idx]
            if b["op"] == group_type:
                group_block(settings, sc, method_idx)
                relaxation_needed = True
            current_block += 1

    # ------------------------------------------------------------------
    # Relax (renumber) the species or types to fill any gaps.  This is
    # only needed when target or block methods were used (reduce does
    # its own complete renumbering internally).  Only applied to file
    # set 0 (the base, non-XANES configuration).
    # ------------------------------------------------------------------
    if relaxation_needed:
        relax_group(settings, group_type, file_set=0)


# ---------------------------------------------------------------------------
# assign_group helper stubs — breadth-first with full documentation,
# to be filled in depth-first.
# ---------------------------------------------------------------------------

def prepare_targets(settings, sc):
    """Convert each target's location to both fractional-abc and direct-xyz.

    Each target on the command line was specified as one of:
      - An atom number (loc_type == 1): use that atom's coordinates.
      - A Cartesian x,y,z point (loc_type == 2): convert to fractional.
      - A fractional a,b,c point (loc_type == 3): convert to Cartesian.

    After this function, every target dict in ``settings.targets`` will
    have ``fract_abc`` (a 1-indexed [None, a, b, c] list) and
    ``direct_xyz`` (a 1-indexed [None, x, y, z] list) populated.

    This data is needed by ``make_min_dist_matrices`` to build the
    atom-to-target distance matrix, and by ``group_target`` to test
    each atom's distance to each target.

    Parameters
    ----------
    settings : ScriptSettings
        Must have ``settings.targets`` populated with loc/loc_type.
    sc : StructureControl
        Provides coordinate arrays and conversion methods.
    """
    print("      Preparing Targets for Inclusion in Distance Matrix")

    for t in settings.targets:
        if t["loc_type"] == 1:
            # Target is an atom number — copy that atom's coordinates.
            atom_idx = t["loc"]
            t["fract_abc"] = list(sc.fract_abc[atom_idx])    # [a, b, c]
            t["direct_xyz"] = list(sc.direct_xyz[atom_idx])  # [x, y, z]

        elif t["loc_type"] == 2:
            # Target given in Cartesian x, y, z.
            xyz = t["loc"]   # [x, y, z] already parsed as floats
            t["direct_xyz"] = list(xyz)

            # Convert xyz → direct_abc (fractional × mag), then ÷ mag → fract.
            # Use direct_xyz2fract_abc which gives pure fractional coords.
            fract = sc.direct_xyz2fract_abc(xyz)
            t["fract_abc"] = fract

        elif t["loc_type"] == 3:
            # Target given in fractional a, b, c.
            # The Perl code divides the input values by mag to get fractional
            # coordinates, meaning the input is actually in direct-space abc
            # (Angstroms projected along lattice vectors).
            abc_input = t["loc"]   # [a, b, c] as floats
            fract = [abc_input[i] / sc.mag[i + 1] for i in range(3)]
            t["fract_abc"] = fract

            # Convert fractional → Cartesian xyz.
            xyz = sc.fract_abc2direct_xyz(fract)
            t["direct_xyz"] = xyz

    settings._targets_prepared = True


def make_min_dist_matrices(settings, sc):
    """Build the minimum-distance matrix between all atoms (and targets).

    Wraps ``sc.create_min_dist_matrix`` using the atom fractional and
    Cartesian coordinates plus any target points.  After this call the
    distance matrix is accessible as ``sc.min_dist`` (atom-atom) and
    the extra rows/columns correspond to target points.

    The minimum-distance matrix accounts for periodic boundary conditions:
    for each pair, it finds the shortest distance considering all lattice
    translation images within the interaction limit distance.

    Also builds ``sc.self_min_dist`` — the distance from each atom to its
    own nearest periodic image.

    Parameters
    ----------
    settings : ScriptSettings
        Provides target coordinates (``targets[i]["fract_abc"]`` etc.)
        and ``num_atoms``.
    sc : StructureControl
        Provides ``fract_abc``, ``direct_xyz``, and the matrix builder.
    """
    # Build 1-indexed target coordinate arrays for the second item group.
    num_targets = len(settings.targets)
    target_fract_abc = [None]
    target_direct_xyz = [None]
    for t in settings.targets:
        target_fract_abc.append(t["fract_abc"])    # [a, b, c] 0-indexed entry
        target_direct_xyz.append(t["direct_xyz"])  # [x, y, z] 0-indexed entry

    # Call StructureControl's matrix builder.  Group 1 = atoms (defaults),
    # group 2 = target points.
    sc.create_min_dist_matrix(
        num_items2=num_targets,
        items2_abc=target_fract_abc,
        items2_xyz=target_direct_xyz,
    )

    settings._min_dist_made = True


def init_track_flag(settings, relation, op):
    """Initialise the tracking array for "alike" grouping operations.

    When the relation is "alike", atoms that match are all given the
    *same* new species/type number.  We need a flag array to record
    which elements (for species) or which element/species pairs (for
    types) have been touched, so that their counts can be incremented
    exactly once after all atoms have been processed.

    Returns
    -------
    track_flag : list
        For species ops: 1-indexed by element, initialised to 0.
        For type ops: 2D, 1-indexed [element][species], initialised to 0.
        Returns an empty list if relation is not "alike".
    """
    if relation != "alike":
        return []

    if op == "species":
        # 1-indexed by element: track_flag[element] = 0 or 1.
        track_flag = [None]
        for element in range(1, settings.num_elements + 1):
            track_flag.append(0)
        return track_flag

    elif op == "types":
        # 2D, 1-indexed: track_flag[element][species] = 0 or 1.
        track_flag = [None]
        for element in range(1, settings.num_elements + 1):
            species_flags = [None]
            for species in range(1, settings.num_species[element] + 1):
                species_flags.append(0)
            track_flag.append(species_flags)
        return track_flag

    return []


def compare_status_and_request(settings, status, zone, relation, op,
                               atom, track_flag, file_set=0):
    """Decide whether to extend the group for a given atom.

    Combines the atom's spatial status (inside=1 or outside=0 the
    target/block zone) with the requested zone (in/out) and relation
    (alike/diff) to decide whether ``extend_group`` should be called.

    The logic is straightforward:
      - If status == 1 (atom is inside) and zone == "in", extend.
      - If status == 0 (atom is outside) and zone == "out", extend.
      - Otherwise, do nothing for this atom.

    Parameters
    ----------
    settings : ScriptSettings
    status : int
        1 if the atom is inside the zone, 0 if outside.
    zone : str
        "in" or "out".
    relation : str
        "alike" or "diff".
    op : str
        "species" or "types".
    atom : int
        1-based atom index.
    track_flag : list
        The tracking array (see ``init_track_flag``).
    file_set : int
        Which XANES file set (0 = base).
    """
    # If the atom is inside the zone and we want "in", or outside and
    # we want "out", then extend the group for this atom.
    if (status == 1 and zone == "in") or (status == 0 and zone == "out"):
        extend_group(settings, op, atom, relation, file_set, track_flag)


def update_from_track_flag(settings, relation, op, track_flag, file_set=0):
    """Increment species/type counts for elements touched by "alike" ops.

    After all atoms have been processed by ``extend_group`` with the
    "alike" relation, the tracking flag array records which elements
    (or element/species pairs) received new alike assignments.  This
    function increments ``num_species`` or ``num_types`` by 1 for each
    flagged entry.

    Parameters
    ----------
    settings : ScriptSettings
    relation : str
        Only acts when "alike".
    op : str
        "species" or "types".
    track_flag : list
        The tracking array from ``init_track_flag``.
    file_set : int
        Which XANES file set (0 = base).
    """
    if relation != "alike":
        return

    if op == "species":
        for element in range(1, settings.num_elements + 1):
            if track_flag[element] == 1:
                settings.num_species[element] += 1
    elif op == "types":
        for element in range(1, settings.num_elements + 1):
            for species in range(1, settings.num_species[element] + 1):
                if track_flag[element][species] == 1:
                    settings.num_types[file_set][element][species] += 1


def extend_group(settings, op, atom, relation, file_set, track_flag):
    """Assign a new species or type number to a single atom.

    This is the core assignment operation called by ``group_target`` and
    ``group_block`` for each atom that passes the spatial test.

    For **species** operations:
      - The atom's species ID is set to (current num_species + 1).
      - If the relation is "diff", ``num_species`` is immediately
        incremented.
      - If the relation is "alike", the track_flag is set so that
        ``num_species`` will be incremented later (once, after all atoms
        sharing this "alike" assignment have been processed).

    For **types** operations:
      - Same logic but for type IDs within the atom's element/species.

    Parameters
    ----------
    settings : ScriptSettings
    op : str
        "species" or "types".
    atom : int
        1-based atom index.
    relation : str
        "alike" or "diff".
    file_set : int
        Which XANES file set (0 = base).
    track_flag : list
        The tracking array from ``init_track_flag``.
    """
    if op == "species":
        # Get this atom's element ID.
        elem_id = settings.atom_element_id[atom]

        # Assign species ID = current count + 1.
        settings.atom_species_id[atom] = settings.num_species[elem_id] + 1

        if relation == "diff":
            # Each "diff" atom gets its own unique species number.
            settings.num_species[elem_id] += 1
        else:
            # "alike" — flag for later bulk increment.
            track_flag[elem_id] = 1

        # Record that this new species has 1 type.  We may need to extend
        # the num_types list for this element if the species index is new.
        elem_types = settings.num_types[file_set][elem_id]
        new_species_idx = settings.num_species[elem_id]
        while len(elem_types) <= new_species_idx:
            elem_types.append(0)
        elem_types[new_species_idx] = 1

    elif op == "types":
        # Get this atom's element and species IDs.
        elem_id = settings.atom_element_id[atom]
        spec_id = settings.atom_species_id[atom]

        # Assign type ID = current count + 1.
        settings.atom_type_id[file_set][atom] = (
            settings.num_types[file_set][elem_id][spec_id] + 1)

        if relation == "diff":
            # Each "diff" atom gets its own unique type number.
            settings.num_types[file_set][elem_id][spec_id] += 1
        else:
            # "alike" — flag for later bulk increment.
            track_flag[elem_id][spec_id] = 1


def group_block(settings, sc, block_idx):
    """Group atoms by a 3D slab (block) in lattice coordinates.

    Tests each atom to see if its direct-abc position falls within the
    block borders.  The block is defined by from/to values along each
    lattice direction (a, b, c).  The special values 'a', 'b', 'c' for
    the "to" border are replaced by the corresponding lattice magnitude.

    Atoms that are inside (or outside, per the zone setting) the block
    are then processed by ``extend_group`` according to the relation
    (alike/diff) and operand (species/types).

    Parameters
    ----------
    settings : ScriptSettings
    sc : StructureControl
    block_idx : int
        Index into ``settings.blocks``.
    """
    block = settings.blocks[block_idx]
    op       = block["op"]
    zone     = block["zone"]
    relation = block["relation"]
    borders  = block["borders"]   # [[from_a,to_a],[from_b,to_b],[from_c,to_c]]

    epsilon = 1e-5

    # Initialise the tracking array for "alike" grouping.
    track_flag = init_track_flag(settings, relation, op)

    # Resolve any symbolic "to" border values.  The special tokens "a",
    # "b", "c" mean "extend to the full lattice magnitude in that
    # direction".  This lets the user say, for example,
    #   -block -abc 0.0 a 0.0 b 0.0 2.5
    # to create a slab that spans the full a and b dimensions but only
    # the first 2.5 Å in c.
    mag_labels = ["a", "b", "c"]
    for dim in range(3):
        if str(borders[dim][1]).lower() == mag_labels[dim]:
            borders[dim][1] = sc.mag[dim + 1]
        else:
            borders[dim][1] = float(borders[dim][1])
        borders[dim][0] = float(borders[dim][0])

    # Check each atom to determine its position relative to the block
    # borders.  direct_abc[atom] is [a, b, c] (0-indexed within entry).
    for atom in range(1, settings.num_atoms + 1):
        abc = sc.direct_abc[atom]   # [a, b, c]

        # An atom is inside the block if its a, b, c coordinates all
        # fall within the from/to range (inclusive, with epsilon tolerance).
        if (abc[0] >= borders[0][0] - epsilon and
            abc[0] <= borders[0][1] + epsilon and
            abc[1] >= borders[1][0] - epsilon and
            abc[1] <= borders[1][1] + epsilon and
            abc[2] >= borders[2][0] - epsilon and
            abc[2] <= borders[2][1] + epsilon):
            status = 1
        else:
            status = 0

        compare_status_and_request(settings, status, zone, relation, op,
                                   atom, track_flag)

    # Finalise "alike" counts.
    update_from_track_flag(settings, relation, op, track_flag)


def group_target(settings, sc, target_idx):
    """Group atoms by proximity to a target point.

    Uses the minimum-distance matrix to find each atom's distance to the
    target point.  Atoms within (or outside, per the zone setting) the
    target radius are processed by ``extend_group`` according to the
    relation (alike/diff) and operand (species/types).

    The target's distance to each atom is stored in the extra rows of the
    minimum-distance matrix (row = num_atoms + target_index).

    Parameters
    ----------
    settings : ScriptSettings
    sc : StructureControl
    target_idx : int
        Index into ``settings.targets``.
    """
    target = settings.targets[target_idx]
    op       = target["op"]
    zone     = target["zone"]
    relation = target["relation"]
    radius   = target["radius"]

    # Initialise the tracking array for "alike" grouping.
    track_flag = init_track_flag(settings, relation, op)

    # Check each atom's minimum distance to the target point.  The
    # target distances live in the extra rows of the min_dist matrix,
    # at row index (num_atoms + 1-based_target_index).  The target_idx
    # here is a 0-based index into settings.targets, so the 1-based
    # matrix row is num_atoms + target_idx + 1.
    target_row = settings.num_atoms + target_idx + 1

    for atom in range(1, settings.num_atoms + 1):
        current_distance = sc.min_dist[target_row][atom]

        # An atom is "inside" the target sphere if its distance is
        # within the specified radius.
        if current_distance <= radius:
            status = 1
        else:
            status = 0

        compare_status_and_request(settings, status, zone, relation, op,
                                   atom, track_flag)

    # Finalise "alike" counts.
    update_from_track_flag(settings, relation, op, track_flag)


def group_reduce(settings, sc, reduce_idx):
    """Group atoms by local-environment similarity (reduce scheme).

    This is the most complex grouping method.  It classifies atoms of the
    same element into species based on how similar their local chemical
    environments are.  The idea is that two atoms which "see" the same
    arrangement of neighbors out to some distance should be treated as
    equivalent (same species) for the purposes of the calculation.

    Overview of the algorithm
    ~~~~~~~~~~~~~~~~~~~~~~~~~
    The reduce scheme works in two major phases:

    **Phase 1 — Build shell configurations.**

    For every atom in the system, we construct a set of concentric
    spherical shells ("levels") describing its local neighborhood:

    1. Start with the full atom list.  Mark the atom itself as excluded
       (level = -1) and all others as unassigned (level = 0).
    2. For level 1, find the *closest* unassigned atom.  Record that
       distance as the "level distance".  Then mark every unassigned atom
       whose distance falls in the range [level_distance,
       level_distance + thick] *and* is within the overall cutoff as
       belonging to level 1.
    3. Repeat for level 2, 3, … up to the requested number of levels,
       each time only considering atoms not yet assigned to a level.

    The result for each atom is:
    - ``level_distance[atom][level]`` — the distance to the nearest
      neighbor that starts that shell.
    - ``reduced_atoms[atom][level]`` — the list of atoms in that shell.

    **Phase 2 — Compare and assign species.**

    Walk the atom list in order.  For each atom that has not yet been
    assigned a species (a "reduction atom", RA):

    1. Give it a new species ID (increment the species counter for its
       element).
    2. Compare it against every subsequent atom (the "comparison atom",
       CA) of the same element.  The CA is considered equivalent if it
       passes **all** of the following tests at every level:

       a. **Distance test** — the absolute difference in level distances
          between RA and CA must not exceed ``tolerance * RA_distance``.
       b. **Count test** — RA and CA must have the same number of
          neighbors in each shell.
       c. **Composition test** — the set of (element, species) pairs of
          the neighbors must be a permutation of each other.  This is
          checked by a greedy matching algorithm: for each RA neighbor,
          search the CA neighbor list for a match; if found, remove it
          from the CA list and continue.  If any RA neighbor cannot be
          matched, the atoms are not equivalent.

    3. If CA passes all tests, assign it the same species as RA.

    After all atoms are processed, the temporary species assignments are
    copied into the permanent arrays and the type counts are re-initialized
    to 1 (since reduce only operates on species, not types).

    Why "reduce"?
    ~~~~~~~~~~~~~
    The name comes from the idea of *reducing* a large number of nominally
    distinct atoms to a smaller set of equivalent species.  In a perfect
    crystal, all atoms of the same element and Wyckoff position would be
    reduced to a single species.  In an amorphous system, the reduction
    may produce many species if local environments differ significantly,
    or very few if the tolerance is loose.

    Diagnostic output
    ~~~~~~~~~~~~~~~~~
    A ``reduceSummary`` file is written with details of each atom's level
    distances and the element names of its shell neighbors.  This is
    invaluable for understanding why two atoms were or were not grouped
    together, and for tuning the reduce parameters.

    Parameters
    ----------
    settings : ScriptSettings
        The global settings object.  Key attributes used:

        - ``reduces[reduce_idx]`` : dict with keys:
            - ``"level"``     : int   — number of shells to build.
            - ``"thick"``     : float — thickness of each shell (Å).
            - ``"cutoff"``    : float — maximum distance to consider (Å).
            - ``"tolerance"`` : float — fractional tolerance for distance
              comparison (e.g. 0.05 means 5%).
            - ``"op"``        : str   — must be ``"species"`` (types not
              supported for reduce).
        - ``num_atoms``, ``num_elements`` : int
        - ``atom_element_id[atom]``   : int (1-indexed)
        - ``atom_species_id[atom]``   : int (1-indexed) — read for NN
          species identification, then overwritten with new assignments.
        - ``atom_element_name[atom]`` : str — used for diagnostic output.
        - ``num_species[element]``    : int (1-indexed) — overwritten.
        - ``num_types[file_set][element][species]`` : reset to 1.

    sc : StructureControl
        Provides ``sc.min_dist[atom1][atom2]`` — the minimum-image
        distance matrix (must already be computed via
        ``make_min_dist_matrices``).

    reduce_idx : int
        Index into ``settings.reduces`` identifying which reduce
        specification to apply.
    """

    r = settings.reduces[reduce_idx]
    num_levels = r["level"]
    thick      = r["thick"]
    cutoff     = r["cutoff"]
    tolerance  = r["tolerance"]
    op         = r["op"]

    # The reduce method currently only supports species grouping.
    if op != "species":
        print("Only can reduce species now.  Aborting")
        sys.exit(1)

    # Open the diagnostic summary file.
    reduce_fh = open("reduceSummary", "w")

    num_atoms    = settings.num_atoms
    num_elements = settings.num_elements

    # ==================================================================
    # Phase 1: Build the shell configuration for every atom.
    #
    # For each atom we determine which other atoms fall in each
    # concentric spherical shell ("level").  The shells are built
    # outward: at each level we find the closest still-unassigned atom,
    # then sweep up everything within [that distance, that distance +
    # thick] that is also within the overall cutoff.
    # ==================================================================

    # reduced_atoms[atom][level] = list of atom indices in that shell.
    # level_distance[atom][level] = distance to the nearest atom that
    #   starts that shell.
    # Both are 1-indexed on atom and level.
    reduced_atoms  = [None] * (num_atoms + 1)
    level_distance = [None] * (num_atoms + 1)

    for atom in range(1, num_atoms + 1):

        # Initialize an array recording which level each other atom is
        # assigned to for the current atom's shells.  0 = unassigned,
        # -1 = the atom itself (excluded from its own shells).
        atom_level = [0] * (num_atoms + 1)  # index 0 unused
        atom_level[atom] = -1  # The current atom is excluded.

        # Allocate storage for this atom's level distances and shell
        # member lists.  Index 0 is a placeholder.
        level_distance[atom] = [None] * (num_levels + 1)
        reduced_atoms[atom]  = [None] * (num_levels + 1)

        # Determine the configuration for each level separately.
        for level in range(1, num_levels + 1):

            # Find the closest unassigned atom (atom_level == 0) to the
            # current atom.  This atom's distance defines the start of
            # the shell for this level.
            closest_atom = 0
            for atom2 in range(1, num_atoms + 1):
                if atom_level[atom2] == 0:
                    if closest_atom == 0:
                        closest_atom = atom2
                    elif (sc.min_dist[atom][atom2] <
                          sc.min_dist[atom][closest_atom]):
                        closest_atom = atom2

            # Record the distance to this level for the current atom.
            level_distance[atom][level] = sc.min_dist[atom][closest_atom]

            # Mark all atoms whose distance from the current atom falls
            # in the range [closest_distance, closest_distance + thick]
            # AND is within the overall cutoff.  These atoms belong to
            # this level's shell.
            closest_dist = sc.min_dist[atom][closest_atom]
            for atom2 in range(1, num_atoms + 1):
                d = sc.min_dist[atom][atom2]
                if (d >= closest_dist and
                        d <= closest_dist + thick and
                        d <= cutoff):
                    atom_level[atom2] = level

        # ---------------------------------------------------------------
        # Save the level assignments for this atom by collecting the
        # atoms assigned to each level into per-level lists.
        # ---------------------------------------------------------------

        # num_reduced[level] counts atoms in each shell.
        num_reduced = [0] * (num_levels + 1)

        for level in range(1, num_levels + 1):
            reduced_atoms[atom][level] = [None]  # 1-indexed list

        for atom2 in range(1, num_atoms + 1):
            for level in range(1, num_levels + 1):
                if atom_level[atom2] == level:
                    num_reduced[level] += 1
                    reduced_atoms[atom][level].append(atom2)

    # ==================================================================
    # Phase 2: Compare shell configurations and assign species.
    #
    # At this point the 3D structure reduced_atoms[atom][level] is
    # complete for each requested level of reduction.  Now we compare
    # the data to determine which atoms are sufficiently similar to be
    # grouped together as the same species.
    # ==================================================================

    # Initialize a temporary atom species ID array (0 = not yet assigned).
    temp_atom_species_id = [0] * (num_atoms + 1)

    # Initialize a temporary counter for the number of species per element.
    temp_num_species = [0] * (num_elements + 1)

    # Clear the type counts — they will be re-initialized to 1 after
    # the new species assignments are finalized.
    settings.num_types = []

    for atom in range(1, num_atoms + 1):

        # If this atom already has a species assignment (from being
        # matched to an earlier atom), skip it.
        if temp_atom_species_id[atom] != 0:
            continue

        # This is a new "reduction atom" (RA).  Increment the species
        # counter for its element and assign it the new species ID.
        elem = settings.atom_element_id[atom]
        temp_num_species[elem] += 1
        temp_atom_species_id[atom] = temp_num_species[elem]

        # --- Diagnostic output for this reduction atom. ---
        reduce_fh.write("---------------------------------------\n")
        reduce_fh.write(
            f"{atom} {settings.atom_element_name[atom]}"
            f"{temp_atom_species_id[atom]}\n"
        )

        # Build lists of the nearest-neighbor element and species IDs
        # for the RA at each level (used for the composition test).
        # ra_elements[level][i] and ra_species[level][i] are 1-indexed
        # on i (matching the reduced_atoms indexing).
        ra_elements = [None] * (num_levels + 1)
        ra_species  = [None] * (num_levels + 1)

        for level in range(1, num_levels + 1):

            # Record the distance to this level in the summary file.
            reduce_fh.write(
                f"Distance to level {level} = "
                f"{level_distance[atom][level]}\n"
            )

            # Build an element-name string for the summary file, and
            # record element/species IDs for the composition test.
            level_atoms_str = ""
            ra_elements[level] = [None]  # 1-indexed
            ra_species[level]  = [None]  # 1-indexed

            # Walk the list of atoms in this level's shell (1-indexed).
            for ri in range(1, len(reduced_atoms[atom][level])):
                ra_atom = reduced_atoms[atom][level][ri]
                ra_elements[level].append(settings.atom_element_id[ra_atom])
                ra_species[level].append(settings.atom_species_id[ra_atom])
                level_atoms_str += (
                    f" {settings.atom_element_name[ra_atom]}"
                )

            # Print the list of atoms in this level.
            reduce_fh.write(f"Atoms in level {level} :  {level_atoms_str}\n")

        # ---------------------------------------------------------------
        # Check each subsequent atom to see if it is sufficiently
        # similar to the current RA according to the requested criteria.
        # ---------------------------------------------------------------

        for atom2 in range(atom + 1, num_atoms + 1):

            # --- Test 0: Same element? ---
            # If the comparison atom (CA) is not the same element as the
            # RA, it cannot possibly be the same species.
            if settings.atom_element_id[atom] != settings.atom_element_id[atom2]:
                continue

            # --- Test 1: Level distances within tolerance? ---
            # Check that the CA has distances to each level that are
            # within (tolerance * RA_level_distance) of the RA distances.
            failed = False
            for level in range(1, num_levels + 1):
                level_dist_diff = abs(
                    level_distance[atom][level] -
                    level_distance[atom2][level]
                )
                if level_dist_diff > tolerance * level_distance[atom][level]:
                    failed = True
                    break
            if failed:
                continue

            # --- Test 2: Same number of neighbors at each level? ---
            # The number of reduced atoms in each shell must match.
            failed = False
            for level in range(1, num_levels + 1):
                # len - 1 because the lists are 1-indexed (index 0 = None).
                ra_count = len(reduced_atoms[atom][level]) - 1
                ca_count = len(reduced_atoms[atom2][level]) - 1
                if ra_count != ca_count:
                    failed = True
                    break
            if failed:
                continue

            # --- Test 3: Identical element/species composition? ---
            # Check for identical element/species combinations at each
            # level.  The order of atoms within a level is not required
            # to match — only the multiset of (element, species) pairs
            # needs to be the same.
            #
            # This is tested by a greedy matching algorithm: for each RA
            # neighbor, search the CA neighbor list for a matching
            # (element, species) pair.  If found, remove it from the CA
            # list (by swapping it to the front and popping).  If every
            # RA neighbor finds a match, the compositions are equivalent.
            # If any RA neighbor fails to find a match, the atoms are
            # not equivalent and we skip to the next CA.
            failed = False
            for level in range(1, num_levels + 1):

                # Build the CA element and species lists for this level.
                # These are consumed (shrunk) during matching, so we
                # build fresh copies each time.
                ca_elements = []
                ca_species  = []
                for ri in range(1, len(reduced_atoms[atom2][level])):
                    ca_atom = reduced_atoms[atom2][level][ri]
                    ca_elements.append(settings.atom_element_id[ca_atom])
                    ca_species.append(settings.atom_species_id[ca_atom])

                # Try to match each RA neighbor against the CA list.
                level_failed = False
                for ri in range(1, len(ra_species[level])):
                    ra_sp = ra_species[level][ri]
                    ra_el = ra_elements[level][ri]

                    # Search the CA list for a matching pair.
                    match_found = False
                    for ci in range(len(ca_species)):
                        if ca_species[ci] == ra_sp and ca_elements[ci] == ra_el:
                            # Match found.  Remove from CA list by
                            # swapping with the first element and
                            # popping the front (matching the Perl
                            # swap-and-shift idiom).
                            ca_species[ci]  = ca_species[0]
                            ca_elements[ci] = ca_elements[0]
                            ca_species.pop(0)
                            ca_elements.pop(0)
                            match_found = True
                            break

                    if not match_found:
                        level_failed = True
                        break

                if level_failed:
                    failed = True
                    break

            if failed:
                continue

            # By this point all tests have been passed and this
            # comparison atom (CA) is to be considered the same
            # element/species as the reduction atom (RA).
            temp_atom_species_id[atom2] = temp_atom_species_id[atom]

    # ==================================================================
    # Finalize: copy temporary results into the permanent arrays.
    # ==================================================================

    # Copy the atom species ID numbers.
    for atom in range(1, num_atoms + 1):
        settings.atom_species_id[atom] = temp_atom_species_id[atom]

    # Copy the number of species for each element.
    for element in range(1, num_elements + 1):
        settings.num_species[element] = temp_num_species[element]

    # Re-initialize the number of types for each species to 1 (since
    # reduce only assigns species, every species starts with one type).
    # num_types is indexed as [file_set][element][species], all 1-indexed.
    # We only populate file_set 0 here; XANES file sets are handled later.
    settings.num_types = [[None]]  # file_set 0, element index 0 placeholder
    for element in range(1, num_elements + 1):
        element_types = [None]  # species index 0 placeholder
        for species in range(1, settings.num_species[element] + 1):
            element_types.append(1)
        settings.num_types[0].append(element_types)

    reduce_fh.close()


def relax_group(settings, group_type, file_set=0):
    """Renumber species or types to fill gaps left by grouping operations.

    After target/block grouping, species or type IDs may have gaps (i.e.
    they are no longer contiguous).  For example, suppose a silicon element
    originally has species Si1 through Si9 assigned across its atoms, but
    after a target-based grouping operation only Si1, Si3, Si4, and Si9
    are still in use.  The sequence has holes at Si2, Si5–Si8.  This
    function compresses ("relaxes") the numbering so that all IDs are
    contiguous starting from 1.

    Why this matters
    ~~~~~~~~~~~~~~~~
    Downstream output routines (basis set printing, potential file lookup,
    etc.) iterate over ``1..num_species[element]`` or
    ``1..num_types[file_set][element][species]`` and expect every ID in
    that range to correspond to at least one atom.  Holes would cause
    empty entries in the output files and waste computational resources
    on non-existent groups.

    Algorithm (three passes)
    ~~~~~~~~~~~~~~~~~~~~~~~~
    The same logical pattern is used for both species and types, but
    types require one additional nesting level (per-species within each
    element).

    **Pass 1 — Flag.**  Build a boolean "flagger" array indexed by the
    current ID values.  Walk the atom list and set ``flagger[id] = 1``
    for every ID that is actually assigned to at least one atom.

    **Pass 2 — Map.**  Walk the flagger array in order.  Maintain a
    running counter that increments only when a flagged (in-use) ID is
    encountered.  Record the counter value at every position — this is
    the *new* compressed ID for that old ID.  Positions that were never
    flagged receive the same counter value as the previous position, but
    that is harmless because no atom will ever look up an unflagged ID.

    **Pass 3 — Apply.**  Walk the atom list again and replace each atom's
    old ID with the mapped (compressed) value.  Also update the
    ``num_species`` or ``num_types`` counts (the final counter value from
    Pass 2 gives the new count).

    Worked example (species)
    ~~~~~~~~~~~~~~~~~~~~~~~~
    Given silicon species after various grouping operations::

        In use:  Si1  Si3  Si4  Si9
        Flagger: [_, 1, 0, 1, 1, 0, 0, 0, 0, 1]   (index 0 unused)
        Map:     [_, 1, 1, 2, 3, 3, 3, 3, 3, 4]
        Result:  Si1→Si1, Si3→Si2, Si4→Si3, Si9→Si4

    Note that even IDs that do not really exist (e.g. Si2, Si7) receive a
    mapping number — it will be the same as the preceding entry.  However,
    that value will never be referenced because no atom carries those IDs.
    Consider index 2: species Si2 does not exist, so it could never be
    mapped to Si1.  Or consider Si7: it does not exist, so it would never
    be mapped to Si3.  But Si9 *does* exist and it will be mapped to Si4.

    Assumptions
    ~~~~~~~~~~~
    - If species are being relaxed, types have not yet been assigned
      (every atom still has type 1 for every file set).
    - If types are being relaxed, species have already been completely
      assigned and are themselves already relaxed.

    Parameters
    ----------
    settings : ScriptSettings
        The global settings object.  The following attributes are read
        and modified in place:

        - ``num_atoms`` : int — total number of atoms.
        - ``num_elements`` : int — total number of distinct elements.
        - ``atom_element_id[atom]`` : int — element index for each atom
          (1-indexed).
        - ``atom_species_id[atom]`` : int — species index for each atom
          (1-indexed).  Modified when ``group_type == "species"``.
        - ``atom_type_id[file_set][atom]`` : int — type index for each
          atom (1-indexed).  Reset to 1 when relaxing species; remapped
          when relaxing types.
        - ``num_species[element]`` : int — number of species per element
          (1-indexed).  Updated when ``group_type == "species"``.
        - ``num_types[file_set][element][species]`` : int — number of
          types per species per element per file set (1-indexed).
          Updated when ``group_type == "types"``.

    group_type : str
        ``"species"`` or ``"types"`` — which level of grouping to relax.
    file_set : int, optional
        Which input-file set (XANES index) to operate on.  Defaults to 0
        (the base file set used for non-XANES calculations).
    """

    if group_type == "species":
        _relax_species(settings, file_set)
    elif group_type == "types":
        _relax_types(settings, file_set)


def _relax_species(settings, file_set):
    """Compress species IDs for every element so they are contiguous from 1.

    This is the simpler of the two relax operations because species are
    indexed only by element (one level of nesting).  After compressing,
    every atom's type ID for the given file set is reset to 1, because
    type assignment has not yet occurred when species are being relaxed.

    See ``relax_group`` for the full algorithm description.
    """

    # --- Pass 1: Flag which species IDs are actually in use. ---
    # species_flagger[element][species] = 1 if at least one atom of that
    # element currently carries that species ID, 0 otherwise.
    # Both dimensions are 1-indexed (index 0 is a placeholder).
    species_flagger = [None]  # placeholder for index 0
    for element in range(1, settings.num_elements + 1):
        # Initialize flags for every possible species of this element.
        flags = [None]  # placeholder for index 0
        for species in range(1, settings.num_species[element] + 1):
            flags.append(0)
        species_flagger.append(flags)

    # Walk the atom list and flag every species that is encountered.
    for atom in range(1, settings.num_atoms + 1):
        elem = settings.atom_element_id[atom]
        spec = settings.atom_species_id[atom]
        species_flagger[elem][spec] = 1

    # --- Pass 2: Build the old-to-new species mapping. ---
    # species_map[element][old_species] = new_species.
    # The mapping counter increments only when a flagged species is found.
    # Unflagged (unused) species get the same map value as the previous
    # species, but that value will never be looked up by any atom.
    species_map = [None]  # placeholder for index 0
    for element in range(1, settings.num_elements + 1):
        # Initialize a counter for this element.  It increments each time
        # we encounter a species that actually exists (flag == 1).
        current_species_id = 0
        mapping = [None]  # placeholder for index 0

        for species in range(1, settings.num_species[element] + 1):
            # Since this species exists, increment the running counter.
            if species_flagger[element][species] == 1:
                current_species_id += 1

            # Record the current counter as the new ID for this old ID,
            # whether or not the old ID was actually in use.
            mapping.append(current_species_id)

        species_map.append(mapping)

        # As a side effect of this loop we can simultaneously fix the
        # number of species of this element — the final counter value
        # is exactly the count of distinct species that remain.
        settings.num_species[element] = current_species_id

    # --- Pass 3: Apply the mapping to every atom. ---
    for atom in range(1, settings.num_atoms + 1):
        elem = settings.atom_element_id[atom]
        old_spec = settings.atom_species_id[atom]
        settings.atom_species_id[atom] = species_map[elem][old_spec]

        # Types have not yet been assigned, so reset to 1.
        settings.atom_type_id[file_set][atom] = 1


def _relax_types(settings, file_set):
    """Compress type IDs for every species of every element so they are
    contiguous from 1.

    This procedure is one level more complicated than the species
    relaxation above.  Instead of remapping species for each element,
    we must remap types for each species of each element — an additional
    level of nesting.

    The algorithm is identical in structure (flag → map → apply), but the
    flagger and map arrays are indexed by ``[species][type]`` and are
    rebuilt from scratch for each element (cleared between elements to
    avoid stale data from a previous element leaking in).

    See ``relax_group`` for the full algorithm description.
    """

    for element in range(1, settings.num_elements + 1):

        # --- Pass 1: Flag which type IDs are in use for this element. ---
        # type_flagger[species][type] = 1 if at least one atom of this
        # element and species carries that type ID, 0 otherwise.
        type_flagger = [None]  # placeholder for index 0
        for species in range(1, settings.num_species[element] + 1):
            flags = [None]  # placeholder for index 0
            for t in range(1, settings.num_types[file_set][element][species] + 1):
                flags.append(0)
            type_flagger.append(flags)

        # Walk the atom list and flag every type that belongs to this element.
        for atom in range(1, settings.num_atoms + 1):
            if settings.atom_element_id[atom] == element:
                spec = settings.atom_species_id[atom]
                t = settings.atom_type_id[file_set][atom]
                type_flagger[spec][t] = 1

        # --- Pass 2: Build the old-to-new type mapping for this element. ---
        # This mapping scheme works just like the one for species above,
        # but with an extra nesting level: we build a separate mapping
        # for each species within this element.
        type_map = [None]  # placeholder for index 0
        for species in range(1, settings.num_species[element] + 1):
            # Initialize a counter for this species.  It increments each
            # time we encounter a type that actually exists (flag == 1).
            current_type_id = 0
            mapping = [None]  # placeholder for index 0

            for t in range(1, settings.num_types[file_set][element][species] + 1):
                # Since this type exists, increment the running counter.
                if type_flagger[species][t] == 1:
                    current_type_id += 1

                # Record the current counter as the new ID for this old ID,
                # whether or not the old ID was actually in use.
                mapping.append(current_type_id)

            type_map.append(mapping)

            # As a side effect, fix the number of types for this species
            # of this element — the final counter value is the count of
            # distinct types that remain.
            settings.num_types[file_set][element][species] = current_type_id

        # --- Pass 3: Apply the mapping to every atom of this element. ---
        for atom in range(1, settings.num_atoms + 1):
            if settings.atom_element_id[atom] == element:
                spec = settings.atom_species_id[atom]
                old_type = settings.atom_type_id[file_set][atom]
                settings.atom_type_id[file_set][atom] = type_map[spec][old_type]

        # The type_flagger and type_map are local to this element iteration
        # and will be rebuilt fresh for the next element (matching the Perl
        # "undef @typeFlagger; undef @typeMap" at the end of each element).


def assign_xanes_types(settings, sc):
    """Create per-atom type assignments for XANES calculations.

    In an X-ray Absorption Near-Edge Structure (XANES) calculation, one
    atom is "excited" — its core orbitals are included explicitly in the
    basis set rather than being orthogonalized out.  The atoms near the
    excited atom need to have distinct potential types so that their
    response to the core hole can be modeled accurately.

    How it works
    ~~~~~~~~~~~~
    For each XANES target atom (the "XA"), this function:

    1. **Copies** the base file-set type assignments (file set 0) into a
       new file set dedicated to this XA.  This means each XANES atom
       gets its own complete set of type assignments, independent of all
       others.

    2. **Walks** the atom list and checks each atom's distance from the
       XA.  If an atom is within ``xanes_radius``, it is given its own
       unique type via ``extend_group("types", ..., "diff", ...)``.
       Atoms outside the radius keep whatever type they had from the
       base file set.

    3. **Relaxes** the type numbering for this file set so that the type
       IDs are contiguous (no gaps), via ``relax_group("types", ...)``.

    The result is that every XANES target atom produces a separate set
    of input files where the atoms near the target each have unique
    potentials, while atoms far from the target share types as usual.

    Selecting the XANES atom list
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The list of atoms to excite is determined in one of two ways:

    - **Explicit** (``-atom`` sub-option within ``-xanes``): The user
      specifies exactly which atom numbers to excite.  These are stored
      in ``settings.xanes_atoms`` during command-line parsing.

    - **Automatic** (default): If no ``-atom`` list is given, the script
      automatically selects one representative atom from each species of
      each element.  Specifically, for each (element, species) pair, the
      *last* atom encountered in the atom list with that combination is
      chosen.  This is handled by ``_prepare_xanes_list()``.

    The automatic mode is the most common case for systematic XANES
    studies: it ensures that every chemically distinct site in the
    system is sampled exactly once, without the user having to manually
    identify the atom numbers.

    Distance matrix note
    ~~~~~~~~~~~~~~~~~~~~
    The Perl version had to handle a half-filled (upper-triangular)
    minimum-distance matrix by conditionally swapping atom indices.  The
    Python ``StructureControl.create_min_dist_matrix`` fills the full
    symmetric matrix, so we can access ``sc.min_dist[xa][atom]``
    directly regardless of index ordering.

    Parameters
    ----------
    settings : ScriptSettings
        Key attributes used:

        - ``xanes_atoms`` : list of int — explicit atom numbers to
          excite (empty list triggers auto-selection).
        - ``xanes_radius`` : float — radius (Å) of the sphere around
          each XA within which atoms get unique types.
        - ``num_atoms``, ``num_elements``, ``num_species`` : counts.
        - ``atom_element_id``, ``atom_species_id`` : per-atom IDs.
        - ``atom_type_id[file_set][atom]`` : modified for each XA
          file set.
        - ``num_types[file_set][element][species]`` : modified.

    sc : StructureControl
        Provides ``sc.min_dist[atom1][atom2]`` for distance lookups.
    """

    # When the -xanes option is given, the default behavior is to select
    # one atom of each species and create a separate set of input files
    # for it.  This can be modified by using the -atom option within
    # -xanes, in which case only the explicitly listed atoms are used.
    # When -atom is given, xanes_atoms is populated during CLI parsing.
    # When it is not given, xanes_atoms is left empty (length 0).
    if len(settings.xanes_atoms) == 0:
        _prepare_xanes_list(settings)

    num_xanes_atoms = len(settings.xanes_atoms)

    # For each atom in the XANES atom list, create a separate set of
    # type assignments.  The XANES atom (XA) is placed at the center
    # of a sphere of radius xanes_radius.  All atoms within that sphere
    # are given their own type.  This is applied independently for each
    # XA, producing a separate input-file set for each one.
    for xa_idx in range(num_xanes_atoms):
        # The file set index for this XANES atom.  File set 0 is the
        # base (non-XANES) configuration; XANES file sets are 1-indexed
        # (xa_idx + 1) matching the Perl convention where $xanesAtom
        # runs 1..$numXanesAtoms and is used directly as the file set.
        file_set = xa_idx + 1
        xa_atom = settings.xanes_atoms[xa_idx]

        # --- Initialize the type IDs for this file set by copying the
        # base (file set 0) assignments. ---
        # Extend atom_type_id to accommodate the new file set if needed.
        while len(settings.atom_type_id) <= file_set:
            settings.atom_type_id.append(
                [None] * (settings.num_atoms + 1)
            )
        for atom in range(1, settings.num_atoms + 1):
            settings.atom_type_id[file_set][atom] = \
                settings.atom_type_id[0][atom]

        # --- Initialize the num_types for this file set by copying the
        # base (file set 0) counts. ---
        while len(settings.num_types) <= file_set:
            settings.num_types.append(None)
        fs_types = [None]  # element index 0 placeholder
        for element in range(1, settings.num_elements + 1):
            elem_types = [None]  # species index 0 placeholder
            for species in range(1, settings.num_species[element] + 1):
                elem_types.append(settings.num_types[0][element][species])
            fs_types.append(elem_types)
        settings.num_types[file_set] = fs_types

        # --- Check each atom's distance from the XA.  If it is within
        # xanes_radius, assign it a unique type ("diff"). ---
        for atom in range(1, settings.num_atoms + 1):
            # The Python min_dist matrix is fully symmetric, so we can
            # reference sc.min_dist[xa_atom][atom] directly (unlike the
            # Perl version which had to conditionally swap indices for
            # the half-filled matrix).
            current_distance = sc.min_dist[xa_atom][atom]

            # Is the current atom inside the XANES sphere?
            if current_distance <= settings.xanes_radius:
                # Assign a unique type to this atom.  "diff" means each
                # atom gets its own type number.  track_flag is unused
                # for "diff" relations, so we pass an empty list.
                extend_group(settings, "types", atom, "diff",
                             file_set, [])

        # Allow the types to be relaxed to a contiguous sequence of
        # numbers for this file set.
        relax_group(settings, "types", file_set=file_set)


def _prepare_xanes_list(settings):
    """Auto-select one representative atom per species for XANES excitation.

    When the user does not explicitly specify which atoms to excite (via
    the ``-atom`` sub-option of ``-xanes``), this function builds the
    list automatically.  The strategy is simple: for each (element,
    species) pair, the *last* atom in the atom list that belongs to that
    pair is selected.

    Why "last"?  In practice the choice is arbitrary — any atom of a
    given species is equally representative.  Using the last one found
    is simply the natural result of overwriting as we walk the atom list
    forward.

    The resulting list is stored in ``settings.xanes_atoms`` as a flat,
    0-indexed Python list of 1-based atom numbers, ordered by element
    first, then species within each element.

    Parameters
    ----------
    settings : ScriptSettings
        Uses ``num_atoms``, ``num_elements``, ``num_species``,
        ``atom_element_id``, ``atom_species_id``.  Sets
        ``xanes_atoms``.
    """

    # Go through the list of atoms.  Record any atom number corresponding
    # to a given element and species.  After going through all atoms,
    # the last atom number for each species of each element will be
    # recorded as the XANES atom.
    xanes_temp = [None]  # 1-indexed by element
    for element in range(1, settings.num_elements + 1):
        elem_list = [None]  # 1-indexed by species
        for species in range(1, settings.num_species[element] + 1):
            elem_list.append(0)
        xanes_temp.append(elem_list)

    for atom in range(1, settings.num_atoms + 1):
        elem = settings.atom_element_id[atom]
        spec = settings.atom_species_id[atom]
        xanes_temp[elem][spec] = atom

    # Parse the 2D xanes_temp structure into a flat list of atom numbers,
    # ordered by element then species.
    settings.xanes_atoms = []
    for element in range(1, settings.num_elements + 1):
        for species in range(1, settings.num_species[element] + 1):
            settings.xanes_atoms.append(xanes_temp[element][species])


def initialize_emu(settings, sc):
    """Create EMU (Evolutionary Multi-objective Utility) configuration files.

    The purpose of this function is to prepare the working directory for an
    EMU optimisation run.  EMU uses particle-swarm optimisation (PSO) to
    explore structural candidates, and requires two small configuration
    files — ``emu.conf`` and ``pso.conf`` — that specify the search
    parameters.

    **What this function does, step by step:**

    1. **Wipe the project directory.**  All files produced by prior makeinput
       steps (structure.dat, olcao.dat, kp files, potentials, etc.) are
       deleted.  Only the original skeleton file ``olcao.skl`` is preserved:
       it is temporarily moved one level up, the project directory is cleared,
       and then it is moved back.

       .. warning::

          This is intentionally destructive.  The ``-emu`` flag is meant to
          be used at the very start of a new EMU campaign, before any
          production files have been generated.  If you have results you
          want to keep, back them up **before** running with ``-emu``.

    2. **Write ``emu.conf``** with default EMU search parameters:

       - Number of swarms: how many independent swarms explore the
         configuration space simultaneously.
       - Number of candidates: the population size within each swarm.
       - Number of unique elements: the count of distinct chemical
         elements that are allowed to vary during the search.
       - Defect cutoff: a distance threshold (in Angstroms) below which
         two candidate atoms are considered to overlap, triggering a
         "defect" penalty.  Candidates with overlapping atoms are
         penalised to steer the search toward physically reasonable
         structures.

    3. **Write ``pso.conf``** with default PSO algorithm parameters:

       - Halting criterion: the convergence rule that decides when a
         swarm has found its optimum (1 = default stagnation-based
         stopping).
       - Number of concurrent threads: how many candidates are evaluated
         in parallel.
       - Number of swarms and candidates: duplicated here because the PSO
         engine reads its own configuration independently of EMU.
       - Number of atoms: taken from the skeleton so the PSO engine knows
         the dimensionality of the search space (3 coordinates per atom).

    Parameters
    ----------
    settings : ScriptSettings
        The command-line / rc-file settings object.  Only
        ``settings.num_atoms`` is read.
    sc : StructureControl
        The structure control instance (not modified, but present for API
        consistency with the other workflow functions).

    Side Effects
    ------------
    - Deletes **all** files in the current project directory except
      ``olcao.skl``.
    - Creates ``emu.conf`` and ``pso.conf`` in the project directory.
    """

    # --- Step 1: Wipe the project directory, preserving only olcao.skl ---
    #
    # The Perl original did:
    #   chdir "..";
    #   system("mv $proj_home/olcao.skl emu_temp");
    #   system("rm -r $proj_home/*");
    #   system("mv emu_temp $proj_home/olcao.skl");
    #   chdir "$proj_home";
    #
    # We replicate the same logic using Python's pathlib / shutil so that
    # the behaviour is platform-consistent and easier to audit.  The
    # project directory is simply the current working directory.

    import shutil
    from pathlib import Path

    proj_dir = Path.cwd()
    skl_name = OLCAO_SKL                      # "olcao.skl"
    skl_path = proj_dir / skl_name

    # Temporarily move the skeleton out of harm's way.
    temp_skl = proj_dir.parent / "emu_temp"
    if skl_path.exists():
        shutil.move(str(skl_path), str(temp_skl))

    # Remove everything inside the project directory.  We iterate over
    # the contents rather than deleting the directory itself so that the
    # working directory remains valid.
    for child in proj_dir.iterdir():
        if child.is_dir():
            shutil.rmtree(child)
        else:
            child.unlink()

    # Restore the skeleton.
    if temp_skl.exists():
        shutil.move(str(temp_skl), str(skl_path))

    # --- Step 2: Write emu.conf with default EMU search parameters ---

    with open(EMU_CONF, "w") as f:
        f.write("Number of swarms               = 1\n")
        f.write("number of candidates           = 10\n")
        f.write("Number of unique elements      = 1\n")
        f.write("Defect cutoff                  = 0.900\n")

    # --- Step 3: Write pso.conf with default PSO algorithm parameters ---

    with open(PSO_CONF, "w") as f:
        f.write("Halting criterion              = 1\n")
        f.write("Number of concurrent threads   = 4\n")
        f.write("Number of swarms               = 1\n")
        f.write("Number of candidates           = 10\n")
        f.write(f"Number of atoms                = {settings.num_atoms}\n")


def print_olcao(settings, sc):
    """Write all OLCAO output files (olcao.dat, structure.dat, scfV, kp, etc.).

    This is the master output orchestrator.  It takes the structure that has
    been loaded, grouped, and (optionally) XANES-assigned, and produces every
    input file that the OLCAO Fortran programs need:

    * ``structure.dat`` — lattice vectors, atomic positions, and type IDs
    * ``olcao.dat``     — the main control file for every OLCAO calculation
      phase (SCF, DOS, bond, SYBD, OPTC, PACS, field, LOEN, etc.)
    * ``scfV.dat``      — starting SCF potential in Gaussian form
    * ``kp-scf.dat`` and ``kp-pscf.dat`` — k-point meshes
    * ``memory``        — estimated peak memory for each OLCAO phase
    * ``datSkl.map``    — mapping between sorted (olcao.dat) atom numbers and
      the original (olcao.skl) atom numbers
    * ``olcao.fract-mi`` and ``olcao.cart-mi`` — re-written skeleton files
      with explicit species assignments (non-XANES only)
    * Job submission scripts (PBS / LSF / Slurm / bash)

    For XANES calculations the loop runs once for the base (no core-hole)
    configuration (file_set 0) and then once more for each XANES target atom
    (file_set 1, 2, ...).  Each XANES file set gets its own ``inputs/``
    directory which is then renamed to a descriptive per-atom directory.

    The steps within each iteration are:

    1. Compute cumulative type counts so that per-element/species type IDs
       can be mapped to a single sequential numbering.
    2. Sort atoms by element, then species, then type (three successive stable
       sorts) and record the old-to-new mapping.
    3. Write ``structure.dat`` (cell vectors and sorted atom positions).
    4. Count the number of atoms belonging to each sequential type.
    5. Write ``olcao.dat`` (basis sets, potentials, control parameters).
    6. Copy k-point files into the ``inputs/`` directory.
    7. Compute and write memory estimates for each OLCAO phase.
    8. (Non-XANES only) Write ``olcao.fract-mi`` and ``olcao.cart-mi``.
    9. Optionally produce PDB and/or CIF visualisation files.
    10. For XANES, move ``inputs/`` into a per-target directory and build a
        ``batchList`` file.
    11. Write the job submission script.

    Parameters
    ----------
    settings : ScriptSettings
        Fully initialised settings object (CLI + rc + cell data).
    sc : StructureControl
        The structure control instance with lattice, atom positions, etc.
    """

    from structure_control import PI, BOHR_RAD

    num_xanes_atoms = len(settings.xanes_atoms)

    # Open a file listing XANES subdirectories for batch submission.
    batch_fh = None
    if num_xanes_atoms > 0:
        batch_fh = open(BATCH_LIST, "w")

    # ---- Step 1: Convert lattice and positions to atomic units --------
    # The OLCAO Fortran programs expect Bohr, not Angstroms.
    _convert_a_to_au(settings, sc)

    # ---- Step 2: Generate k-point meshes (once for all file sets) -----
    _make_kp(settings, sc)

    # ---- Step 3: Contract the basis set for each unique element -------
    _contract_basis(settings, sc)

    # ---- Step 4: Copy potential information from the database ----------
    _obtain_pot_info(settings, sc)

    # ---- Step 5: Loop over file sets (base + XANES targets) -----------
    # The project name is the basename of the current working directory,
    # matching the Perl ``$proj_home``.
    proj_home = os.path.basename(os.getcwd())

    for h in range(0, num_xanes_atoms + 1):
        if h == 0:
            print("      Working on standard input files")
        else:
            print(f"      Working on xanes input file set {h}")

        # Create the inputs directory for this file set.
        os.makedirs(INPUTS_DIR, exist_ok=True)

        # 5a. Cumulative type counts.
        cumulative_num_types, total_num_types = _get_cumulative_types(settings, h)

        # 5b. Sort atoms by element → species → type.
        (sorted_elem_id, sorted_spec_id, sorted_type_id,
         sorted_fract_abc, sorted_direct_abc, sorted_direct_xyz,
         sorted_xanes_atoms, index_map) = _sort_atoms(
            settings, sc, h, cumulative_num_types)

        # 5c. Write structure.dat.
        _print_structure(settings, sc, h, sorted_elem_id, sorted_spec_id,
                         sorted_type_id, sorted_direct_xyz,
                         cumulative_num_types, total_num_types)

        # 5d. Count atoms per sequential type.
        num_atoms_of_type = _get_num_atoms_of_type(
            settings, h, sorted_elem_id, sorted_spec_id,
            sorted_type_id, cumulative_num_types, total_num_types)

        # 5e. Write olcao.dat (the main control file).
        (num_total_core_states, num_total_vale_states,
         num_core_xanes, xanes_control,
         max_num_atom_alphas, pot_dim, nuc_charge) = _print_olcao_input(
            settings, sc, h, sorted_elem_id, sorted_spec_id,
            sorted_type_id, sorted_xanes_atoms, cumulative_num_types,
            total_num_types, num_atoms_of_type)

        # 5f. Copy k-point files into the inputs directory.
        _get_kp()

        # 5g. Compute memory estimates for each basis type
        # and write them.  In density mode the k-point count
        # is not known at this stage so we skip memory
        # estimation entirely.
        if not settings.use_kp_density:
            _compute_and_print_mem(
                settings, sc, h,
                num_total_core_states,
                num_total_vale_states,
                total_num_types,
                num_atoms_of_type,
                max_num_atom_alphas,
                pot_dim, nuc_charge)

        # 5h. (Non-XANES only) Write olcao.mi skeleton files.
        if h == 0:
            _make_olcao_mi(settings, sc, sorted_elem_id, sorted_spec_id,
                           sorted_fract_abc, sorted_direct_xyz,
                           OLCAO_FMI)
            _make_olcao_mi(settings, sc, sorted_elem_id, sorted_spec_id,
                           sorted_fract_abc, sorted_direct_xyz,
                           OLCAO_CMI)

        # 5i. Determine the name for this file set's output directory.
        if h == 0:
            current_name = proj_home
        else:
            xa = sorted_xanes_atoms[h]
            current_name = (
                f"{settings.element_list[sorted_elem_id[xa]]}"
                f"{sorted_spec_id[xa]}"
                f"_{settings.xanes_atoms[h - 1]}")

        # 5j. Optionally produce structure visualisation files.
        if settings.pdb == 1:
            _make_pdb(sc, h, current_name, sorted_elem_id,
                      sorted_spec_id, sorted_fract_abc, sorted_direct_xyz)
        if settings.cif == 1:
            _make_cif(sc, settings, h, current_name, sorted_elem_id,
                      sorted_spec_id, sorted_type_id, sorted_fract_abc)

        # 5k. For XANES, move inputs into a per-target directory.
        if settings.xanes == 1:
            current_dir_name = current_name
            if os.path.exists(current_dir_name):
                shutil.rmtree(current_dir_name)
            os.makedirs(current_dir_name, exist_ok=True)
            shutil.move(INPUTS_DIR, current_dir_name)
            if settings.cif == 1:
                cif_file = current_name + CIF_EXT
                if os.path.exists(cif_file):
                    shutil.move(cif_file, current_dir_name)
            if settings.pdb == 1:
                pdb_file = current_name + PDB_EXT
                if os.path.exists(pdb_file):
                    shutil.move(pdb_file, current_dir_name)
            if h > 0 and batch_fh is not None:
                batch_fh.write(f"{current_dir_name}   0\n")

        # 5l. Write the submission script.
        if settings.xanes == 1:
            _make_sub_file(settings, h, current_name, current_dir_name,
                           proj_home, num_core_xanes, xanes_control)
        else:
            _make_sub_file(settings, h, current_name, ".",
                           proj_home, num_core_xanes, xanes_control)

        # Store per-file-set data that print_summary will need.
        # We save it on settings so it persists after print_olcao returns.
        if not hasattr(settings, '_summary_data'):
            settings._summary_data = {}
        num_electrons, num_states_used = _get_elec_and_states(
            settings, nuc_charge, num_total_core_states,
            num_total_vale_states)
        settings._summary_data[h] = {
            'num_total_vale_states': list(num_total_vale_states),
            'num_total_core_states': list(num_total_core_states),
            'pot_dim': pot_dim,
            'num_electrons': num_electrons,
            'num_states_used': list(num_states_used),
            'total_num_types': total_num_types,
            'sorted_xanes_atoms': list(sorted_xanes_atoms),
        }

    # Store the project home name for print_summary.
    settings._proj_home = proj_home

    # Close the batch list if we opened it.
    if batch_fh is not None:
        batch_fh.close()


# ---------------------------------------------------------------------------
# print_olcao helper functions
# ---------------------------------------------------------------------------

def _convert_a_to_au(settings, sc):
    """Convert lattice vectors and atomic positions from Angstroms to Bohr.

    The OLCAO Fortran programs work exclusively in atomic units (Bohr radii).
    The structure control object stores everything in Angstroms after reading
    the olcao.skl skeleton.  This function divides all lengths by the Bohr
    radius (0.5291772180 Å) to convert to atomic units in-place.

    Specifically, the following quantities are converted:
      * ``sc.mag[1..3]``       — lattice vector magnitudes
      * ``sc.real_lattice``    — 3×3 lattice vectors
      * ``sc.direct_abc``      — Cartesian-scaled fractional positions
      * ``sc.direct_xyz``      — Cartesian positions

    Note: ``sc.fract_abc`` (pure fractional coordinates) is dimensionless
    and does NOT need conversion.

    Parameters
    ----------
    settings : ScriptSettings
        Not modified; only ``settings.num_atoms`` is read.
    sc : StructureControl
        Modified in-place.
    """
    from structure_control import BOHR_RAD

    for abc in range(1, 4):
        sc.mag[abc] /= BOHR_RAD
        for xyz in range(1, 4):
            sc.real_lattice[abc][xyz] /= BOHR_RAD

    for atom in range(1, settings.num_atoms + 1):
        for axis in range(1, 4):
            sc.direct_abc[atom][axis] /= BOHR_RAD
            sc.direct_xyz[atom][axis] /= BOHR_RAD


def _make_kp(settings, sc):
    """Generate k-point files for SCF and post-SCF meshes.

    Two modes are supported:

    **Mesh mode** (default, ``use_kp_density`` is False):
      For each k-point group (1 = SCF, 2 = post-SCF):
        1. Writes ``kpSpecs.dat`` for the ``makekpoints`` program.
        2. Runs ``makekpoints`` to produce an explicit list.
        3. Counts the generated k-points (lines minus 9-line header).
        4. Moves the output to ``.inputTemp/kp-scf.dat`` or
           ``.inputTemp/kp-pscf.dat``.

    **Density mode** (``use_kp_density`` is True):
      Bypasses ``makekpoints`` entirely.  For each group, writes
      a style-code-2 k-point file directly into ``.inputTemp/``.
      The file contains the density value and shift; OLCAO
      computes the per-axis mesh at runtime from the reciprocal
      cell geometry.

    Parameters
    ----------
    settings : ScriptSettings
    sc : StructureControl
    """

    kp_group_file = ["", KP_SCF, KP_PSCF]  # 1-indexed
    settings.kp_num = [0, 0, 0]

    if settings.use_kp_density:
        # Density mode: write style-code-2 files directly.
        print("      Generating KPoints (density mode).")

        # Extract the point group operations from the space
        # group database.  These are embedded in the kpoint
        # file so that OLCAO can perform IBZ reduction at
        # runtime without needing the database itself.
        point_ops = _extract_point_ops(settings)

        for kp_group in range(1, 3):
            dest = os.path.join(
                INPUT_TEMP, kp_group_file[kp_group])
            _write_density_kp_file(
                dest,
                settings.kp_density[kp_group],
                settings.kp_shift,
                point_ops,
                settings.kp_intg_code[kp_group])
        # No kpSpecs.dat is produced in density mode.
    else:
        # Mesh mode: call the makeKPoints executable.
        import subprocess
        print("      Generating KPoints.")

        kp_mesh = [
            None,
            settings.kp_mesh_scf,
            settings.kp_mesh_pscf]

        for kp_group in range(1, 3):
            # Write the makekpoints input file.
            _print_kp_in_file(
                settings, sc,
                kp_mesh[kp_group], kp_group)

            # Run the makekpoints executable.
            cmd = (
                f"{settings.makekpoints_exec}"
                f" {settings.print_bz}"
                f" {settings.scale_factor}")
            subprocess.run(cmd, shell=True, check=True)

            # Count k-points (output has a 9-line header).
            with open(KP_OUT_FILE, "r") as f:
                num_lines = sum(1 for _ in f)
            settings.kp_num[kp_group] = num_lines - 9

            # Move output to .inputTemp.
            _safe_move(
                KP_OUT_FILE,
                os.path.join(
                    INPUT_TEMP,
                    kp_group_file[kp_group]))

        # Archive the last kpSpecs.dat.
        _safe_move(
            KP_IN_FILE,
            os.path.join(INPUT_TEMP, KP_IN_FILE))


def _extract_point_ops(settings):
    """Extract point group operations from the space group database.

    The space group database file has this format:
      - Line 1: description (e.g. ``"F Fm3~m"``)
      - Line 2: root space group number, sub-number
      - Line 3: numSpaceOps, numShifts
      - Then for each of the numSpaceOps operations:
          blank line, 3x3 rotation matrix (3 lines),
          translation vector (1 line)

    The first ``numSpaceOps / numShifts`` operations are
    the pure point group operations (no extra translations
    from centering shifts).  We extract only these.

    Parameters
    ----------
    settings : ScriptSettings
        Must have ``space_db`` and ``space_group_name`` set.

    Returns
    -------
    list of list[list[float]]
        Each element is a 3x3 matrix (list of 3 rows, each
        row a list of 3 floats).
    """
    sg_path = os.path.join(
        settings.space_db, settings.space_group_name)
    with open(sg_path, "r") as f:
        # Skip description line.
        f.readline()
        # Skip root space group number line.
        f.readline()
        # Read number of space ops and number of shifts.
        parts = f.readline().split()
        num_space_ops = int(parts[0])
        num_shifts = int(parts[1])
        num_point_ops = num_space_ops // num_shifts

        # Read each point group operation.
        point_ops = []
        for i in range(num_point_ops):
            # Skip blank line between operations.
            f.readline()
            # Read the 3x3 rotation matrix (3 lines).
            matrix = []
            for row in range(3):
                vals = f.readline().split()
                matrix.append(
                    [float(v) for v in vals[:3]])
            # Skip the translation line.
            f.readline()
            point_ops.append(matrix)

    return point_ops


def _write_density_kp_file(dest_path, density,
                           kp_shift, point_ops,
                           intg_code=0):
    """Write a style-code-2 k-point file for density mode.

    This file is read by ``readKPoints`` in ``kpoints.f90``.
    It tells OLCAO to compute the per-axis k-point mesh
    internally using ``computeAxialKPoints``, based on the
    minimum k-point volume density and the reciprocal cell
    geometry.  It also embeds the point group operations so
    that OLCAO can perform IBZ symmetry reduction at runtime.

    The file format matches the style-code-2 branch of
    ``readKPoints``:
      - ``KPOINT_STYLE_CODE`` = 2
      - ``KPOINT_INTG_CODE``  = integration method
        (0 = Gaussian/histogram, 1 = LAT, etc.)
      - ``MIN_KP_LINE_DENSITY`` = the volume density
        (label is historical; the value is kpoints per
        unit reciprocal-space volume in Bohr^-3)
      - ``KP_SHIFT_A_B_C`` = the a, b, c shift values
      - ``NUM_POINT_OPS`` = number of point group operations
      - ``POINT_OPS`` label, then one 3x3 matrix per
        operation (3 lines each, blank-line separated)

    Parameters
    ----------
    dest_path : str
        Full path to the output file (e.g.
        ``.inputTemp/kp-scf.dat``).
    density : float
        The minimum k-point volume density (kpoints per
        unit reciprocal-space volume, Bohr^-3).
    kp_shift : str
        Space-separated shift values
        (e.g. ``"0.5 0.5 0.5"``).
    point_ops : list of list[list[float]]
        Point group rotation matrices from
        ``_extract_point_ops``.
    intg_code : int, optional
        K-point integration method code.  0 = Gaussian
        /histogram (default), 1 = LAT (linear analytic
        tetrahedron).  Higher integers are reserved for
        future methods.
    """
    with open(dest_path, "w") as f:
        f.write("KPOINT_STYLE_CODE\n")
        f.write("2\n")
        f.write("KPOINT_INTG_CODE\n")
        f.write(f"{intg_code}\n")
        f.write("MIN_KP_LINE_DENSITY\n")
        f.write(f"{density}\n")
        f.write("KP_SHIFT_A_B_C\n")
        f.write(f"{kp_shift}\n")
        f.write("NUM_POINT_OPS\n")
        f.write(f"{len(point_ops)}\n")
        f.write("POINT_OPS\n")
        for i, matrix in enumerate(point_ops):
            if i > 0:
                f.write("\n")
            for row in matrix:
                f.write(f"  {row[0]:12.8f}"
                        f"  {row[1]:12.8f}"
                        f"  {row[2]:12.8f}\n")


def _print_kp_in_file(settings, sc, kp_mesh_request, kp_group):
    """Write the ``kpSpecs.dat`` input file for the makekpoints program.

    The file contains:
      1. The real-space lattice vectors (3×3, in Bohr).
      2. The symmetry operations for this space group (read from the
         space group database).
      3. The k-point mesh dimensions (e.g. ``4 4 4``).
      4. The k-point shift string (e.g. ``0.5 0.5 0.5``).
      5. A gamma flag: 1 if this is a single-k-point gamma calculation
         (mesh is 1×1×1 AND the user explicitly requested it), 0 otherwise.

    Parameters
    ----------
    settings : ScriptSettings
    sc : StructureControl
    kp_mesh_request : list of int
        Three integers [nx, ny, nz] for the mesh in this group.
    kp_group : int
        1 = SCF, 2 = PSCF.
    """

    # Determine if this is a gamma-only calculation.  The condition is
    # that the mesh is 1×1×1 AND the user explicitly set the k-points
    # on the command line (captured by kp_note being "(Gamma)").
    if (kp_mesh_request == [1, 1, 1] and
            settings.kp_note[kp_group] == "(Gamma)"):
        do_gamma = 1
    else:
        do_gamma = 0

    with open(KP_IN_FILE, "w") as kp:
        # Lattice vectors.
        for abc in range(1, 4):
            kp.write(f"{sc.real_lattice[abc][1]} "
                     f"{sc.real_lattice[abc][2]} "
                     f"{sc.real_lattice[abc][3]}\n")

        # Space group symmetry operations from the database.
        sg_path = os.path.join(settings.space_db, settings.space_group_name)
        with open(sg_path, "r") as sg:
            kp.write(sg.read())

        # Mesh parameters.
        kp.write(f" {kp_mesh_request[0]} {kp_mesh_request[1]}"
                 f" {kp_mesh_request[2]}\n")
        kp.write(f"{settings.kp_shift}\n")
        kp.write(f"{do_gamma}\n")


def _contract_basis(settings, sc):
    """Contract basis sets for every unique element/species pair.

    For each element/species combination, this function:
      1. Determines the correct basis set file name (possibly modified by
         ``-subbasis`` substitutions on the command line).
      2. Copies the contraction input from the atomic basis database into
         ``.inputTemp/contract.dat`` and runs the external ``contract``
         program.
      3. Renames the output files to include the element name.
      4. If XANES calculations are requested, repeats the contraction for
         the "nocore" variant of the basis set (where core orbitals are
         explicitly included in the valence so that the core hole can be
         created).

    The resulting ``.wf`` (wave function), ``.dat`` (contraction data),
    and ``.out`` (contraction log) files live in ``.inputTemp/`` and are
    referenced later when writing ``olcao.dat``.

    The function also records the file names for each element/species pair
    in ``settings.basis_files`` and ``settings.basis_files_nc``.

    Parameters
    ----------
    settings : ScriptSettings
    sc : StructureControl
    """
    import subprocess

    print("      Contracting Elemental Wave Functions")

    num_xanes_atoms = len(settings.xanes_atoms)

    # Parse basis substitutions: each entry may be "si" (element only) or
    # "si2" (element + species number).
    basis_sub_element = []
    basis_sub_species = []
    for idx, out_val in enumerate(settings.basis_sub_out):
        import re
        m = re.match(r'^([a-zA-Z]+)(\d+)?$', out_val)
        if m:
            basis_sub_element.append(m.group(1))
            basis_sub_species.append(int(m.group(2)) if m.group(2) else 0)
        else:
            basis_sub_element.append(out_val)
            basis_sub_species.append(0)

    # Enter .inputTemp for easier writing.
    orig_dir = os.getcwd()
    os.chdir(INPUT_TEMP)

    # Initialise the basis file name arrays (1-indexed [element][species]).
    settings.basis_files = [None]
    settings.basis_files_nc = [None]
    for element in range(1, settings.num_elements + 1):
        elem_files = [None]
        elem_files_nc = [None]
        for species in range(1, settings.num_species[element] + 1):
            # Default basis file names.
            elem_name = settings.element_list[element]
            curr_basis = f"contract1_{elem_name}"
            curr_basis_nc = f"nocore_contract1_{elem_name}"

            # Apply any substitutions.
            for s_idx in range(len(settings.basis_sub_out)):
                if (basis_sub_element[s_idx].lower() == elem_name.lower() and
                        (basis_sub_species[s_idx] == species or
                         basis_sub_species[s_idx] == 0)):
                    in_num = settings.basis_sub_in[s_idx]
                    curr_basis = f"contract{in_num}_{elem_name}"
                    curr_basis_nc = f"nocore_contract{in_num}_{elem_name}"

            # Check if this basis file has already been contracted.
            if not os.path.exists(f"{curr_basis}.wf"):
                # Copy the contraction input file.
                db_dir = os.path.join(settings.atomic_bdb, elem_name)
                if settings.no_core == 0:
                    shutil.copy2(os.path.join(db_dir, curr_basis),
                                 "contract.dat")
                else:
                    shutil.copy2(os.path.join(db_dir, curr_basis_nc),
                                 "contract.dat")

                # Copy the isolated-atom potential function.
                shutil.copy2(os.path.join(db_dir, "gauss.fit"), "gauss.fit")

                # Run the contraction program.
                cmd = f"{settings.contract_exec} {settings.do_basis_vis}"
                subprocess.run(cmd, shell=True, check=True)

                # Rename the results files.
                _safe_move("contract.out", f"{curr_basis}.out")
                _safe_move("contract.dat", f"{curr_basis}.dat")
                _safe_move("waveFn.dat", f"{curr_basis}.wf")
                if settings.do_basis_vis == 1:
                    _safe_move("waveFn.plot", f"{curr_basis}.plot")

                # Repeat for the nocore variant if XANES is active.
                if num_xanes_atoms > 0:
                    shutil.copy2(os.path.join(db_dir, curr_basis_nc),
                                 "contract.dat")
                    subprocess.run(
                        f"{settings.contract_exec} 0",
                        shell=True, check=True)
                    _safe_move("contract.out", f"{curr_basis_nc}.out")
                    _safe_move("contract.dat", f"{curr_basis_nc}.dat")
                    _safe_move("waveFn.dat", f"{curr_basis_nc}.wf")

            # Record the file names for this element/species pair.
            elem_files.append(f"{curr_basis}.wf")
            elem_files_nc.append(f"{curr_basis_nc}.wf")

        settings.basis_files.append(elem_files)
        settings.basis_files_nc.append(elem_files_nc)

    # Clean up the temporary potential file.
    if os.path.exists("gauss.fit"):
        os.remove("gauss.fit")

    # Return to the project directory.
    os.chdir(orig_dir)


def _obtain_pot_info(settings, sc):
    """Copy potential and coefficient files from the atomic potential database.

    For each element/species pair, the appropriate ``pot<N>_<elem>`` and
    ``coeff<N>_<elem>`` files are copied from ``$OLCAO_DATA/atomicPDB/<elem>/``
    into ``.inputTemp/``.  Substitutions from the ``-subpot`` command-line
    option are honoured (same element/species matching logic as for basis
    substitutions).

    The resulting file names for each element/species are recorded in
    ``settings.pot_files`` and ``settings.coeff_files``.

    Parameters
    ----------
    settings : ScriptSettings
    sc : StructureControl
    """
    import re

    # Parse potential substitutions.
    pot_sub_element = []
    pot_sub_species = []
    for out_val in settings.pot_sub_out:
        m = re.match(r'^([a-zA-Z]+)(\d+)?$', out_val)
        if m:
            pot_sub_element.append(m.group(1))
            pot_sub_species.append(int(m.group(2)) if m.group(2) else 0)
        else:
            pot_sub_element.append(out_val)
            pot_sub_species.append(0)

    settings.pot_files = [None]
    settings.coeff_files = [None]
    for element in range(1, settings.num_elements + 1):
        elem_pots = [None]
        elem_coeffs = [None]
        elem_name = settings.element_list[element]
        for species in range(1, settings.num_species[element] + 1):
            curr_pot = "pot1"
            curr_coeff = "coeff1"

            for s_idx in range(len(settings.pot_sub_out)):
                if (pot_sub_element[s_idx].lower() == elem_name.lower() and
                        (pot_sub_species[s_idx] == species or
                         pot_sub_species[s_idx] == 0)):
                    in_num = settings.pot_sub_in[s_idx]
                    curr_pot = f"pot{in_num}"
                    curr_coeff = f"coeff{in_num}"

            # Copy files if they haven't been copied yet.
            pot_tagged = f"{curr_pot}_{elem_name}"
            coeff_tagged = f"{curr_coeff}_{elem_name}"
            pot_dst = os.path.join(INPUT_TEMP, pot_tagged)
            if not os.path.exists(pot_dst):
                db_dir = os.path.join(settings.atomic_pdb, elem_name)
                shutil.copy2(os.path.join(db_dir, curr_pot), pot_dst)
                shutil.copy2(os.path.join(db_dir, curr_coeff),
                             os.path.join(INPUT_TEMP, coeff_tagged))

            elem_pots.append(pot_tagged)
            elem_coeffs.append(coeff_tagged)

        settings.pot_files.append(elem_pots)
        settings.coeff_files.append(elem_coeffs)


def _get_cumulative_types(settings, file_set):
    """Compute cumulative type counts through elements and species.

    The OLCAO programs use a single sequential numbering for types across
    all elements and species.  For example, if element 1 has 3 types and
    element 2 has 2 types, the sequential type IDs are 1-3 for element 1
    and 4-5 for element 2.  The "cumulative" count for an element/species
    pair is the number of types that precede it in this sequence.

    Parameters
    ----------
    settings : ScriptSettings
    file_set : int
        Which file set (0 = base, 1+ = XANES target).

    Returns
    -------
    cumulative_num_types : list of lists
        ``cumulative_num_types[element][species]`` = number of types
        preceding this element/species pair.  1-indexed.
    total_num_types : int
        Total number of types across all elements/species for this file set.
    """
    counter = 0
    total = 0
    cumulative = [None]
    for element in range(1, settings.num_elements + 1):
        elem_cum = [None]
        for species in range(1, settings.num_species[element] + 1):
            elem_cum.append(counter)
            nt = settings.num_types[file_set][element][species]
            counter += nt
            total += nt
        cumulative.append(elem_cum)
    return cumulative, total


def _sort_atoms(settings, sc, file_set, cumulative_num_types):
    """Sort atoms by element, species, and type; produce the index map.

    Atoms must be ordered by element → species → type for the OLCAO
    Fortran programs.  This function performs three successive stable sorts
    (by type, then species, then element) and produces sorted versions of
    all key per-atom arrays plus the old↔new index mapping.

    The Perl version used StructureControl::stableSort with in-place swap
    tracking.  Here we use Python's built-in ``sorted()`` which is already
    stable (Timsort), composing the three sorts with a compound key.

    Parameters
    ----------
    settings : ScriptSettings
    sc : StructureControl
    file_set : int
    cumulative_num_types : list of lists

    Returns
    -------
    sorted_elem_id, sorted_spec_id, sorted_type_id : lists
        1-indexed (index 0 = 0 sentinel).
    sorted_fract_abc, sorted_direct_abc, sorted_direct_xyz : lists
        1-indexed, each entry is a 4-element list [None, x, y, z].
    sorted_xanes_atoms : list
        Maps file_set → sorted atom number for the XANES target.
    index_map : list
        ``index_map[sorted_pos] = original_atom_number``.
    """
    num_atoms = settings.num_atoms

    # Build a list of (original_atom, sort_key) and sort it.
    # Sort key: (element_id, species_id, type_id) — primary, secondary, tertiary.
    atoms = list(range(1, num_atoms + 1))
    atoms.sort(key=lambda a: (
        settings.atom_element_id[a],
        settings.atom_species_id[a],
        settings.atom_type_id[file_set][a]
    ))

    # Build the index map: index_map[new_pos] = old_atom_number.
    # Position 0 is a sentinel.
    index_map = [0]  # index 0
    for a in atoms:
        index_map.append(a)

    # Build sorted arrays.
    sorted_elem_id = [0]
    sorted_spec_id = [0]
    sorted_type_id = [0]
    sorted_fract_abc = [[None, 0.0, 0.0, 0.0]]
    sorted_direct_abc = [[None, 0.0, 0.0, 0.0]]
    sorted_direct_xyz = [[None, 0.0, 0.0, 0.0]]

    for new_pos in range(1, num_atoms + 1):
        old = index_map[new_pos]
        sorted_elem_id.append(settings.atom_element_id[old])
        sorted_spec_id.append(settings.atom_species_id[old])
        sorted_type_id.append(settings.atom_type_id[file_set][old])
        sorted_fract_abc.append([
            None,
            sc.fract_abc[old][1],
            sc.fract_abc[old][2],
            sc.fract_abc[old][3],
        ])
        sorted_direct_abc.append([
            None,
            sc.direct_abc[old][1],
            sc.direct_abc[old][2],
            sc.direct_abc[old][3],
        ])
        sorted_direct_xyz.append([
            None,
            sc.direct_xyz[old][1],
            sc.direct_xyz[old][2],
            sc.direct_xyz[old][3],
        ])

    # Map original XANES atom numbers to their sorted positions.
    num_xanes = len(settings.xanes_atoms)
    sorted_xanes_atoms = [0] * (num_xanes + 1)
    if file_set == 0:
        pass  # all zeros
    else:
        target = settings.xanes_atoms[file_set - 1]  # 0-indexed list → 1-based atom
        for new_pos in range(1, num_atoms + 1):
            if index_map[new_pos] == target:
                sorted_xanes_atoms[file_set] = new_pos
                break

    # Write the datSkl.map file.
    map_path = os.path.join(INPUTS_DIR, DAT_SKL_MAP)
    with open(map_path, "w") as f:
        f.write("       DAT#     SKELETON#\n")
        for new_pos in range(1, num_atoms + 1):
            f.write(f"{new_pos:10d} {index_map[new_pos]:10d}\n")

    return (sorted_elem_id, sorted_spec_id, sorted_type_id,
            sorted_fract_abc, sorted_direct_abc, sorted_direct_xyz,
            sorted_xanes_atoms, index_map)


def _get_num_atoms_of_type(settings, file_set, sorted_elem_id,
                           sorted_spec_id, sorted_type_id,
                           cumulative_num_types, total_num_types):
    """Count the number of atoms belonging to each sequential type.

    Parameters
    ----------
    settings : ScriptSettings
    file_set : int
    sorted_elem_id, sorted_spec_id, sorted_type_id : lists (1-indexed)
    cumulative_num_types : list of lists (1-indexed)
    total_num_types : int

    Returns
    -------
    num_atoms_of_type : list
        ``num_atoms_of_type[seq_type]`` = count.  1-indexed.
    """
    nat = [0] * (total_num_types + 1)  # 1-indexed

    for atom in range(1, settings.num_atoms + 1):
        seq_type = (sorted_type_id[atom] +
                    cumulative_num_types[sorted_elem_id[atom]]
                                        [sorted_spec_id[atom]])
        nat[seq_type] += 1

    return nat


def _print_structure(settings, sc, file_set, sorted_elem_id, sorted_spec_id,
                     sorted_type_id, sorted_direct_xyz,
                     cumulative_num_types, total_num_types):
    """Write the ``structure.dat`` file.

    This file tells the OLCAO Fortran programs where the atoms are and
    how they are typed.  It contains:

    * ``CELL_VECTORS`` — the 3×3 real-space lattice matrix (in Bohr).
    * ``NUM_ATOM_SITES`` — the total number of atoms.
    * ``NUM_TYPE_X_Y_Z_ELEM`` — for each atom: sequential number,
      sequential type ID, Cartesian x/y/z (Bohr), and element symbol.
    * ``NUM_POTENTIAL_SITES`` — same block repeated for potential sites
      (in OLCAO the atom sites and potential sites are identical).

    Parameters
    ----------
    settings : ScriptSettings
    sc : StructureControl
    file_set : int
    sorted_elem_id, sorted_spec_id, sorted_type_id : lists (1-indexed)
    sorted_direct_xyz : list of lists (1-indexed)
    cumulative_num_types : list of lists (1-indexed)
    total_num_types : int
    """
    path = os.path.join(INPUTS_DIR, STRUCTURE)
    with open(path, "w") as f:
        # Cell vectors.
        f.write("CELL_VECTORS\n")
        for abc in range(1, 4):
            f.write(f"{sc.real_lattice[abc][1]:26.18f}"
                    f"{sc.real_lattice[abc][2]:26.18f}"
                    f"{sc.real_lattice[abc][3]:26.18f}\n")

        # Number of atom sites.
        f.write("NUM_ATOM_SITES\n")
        f.write(f"{settings.num_atoms:<5d}\n")

        # Atom positions.
        f.write("NUM_TYPE_X_Y_Z_ELEM\n")
        for atom in range(1, settings.num_atoms + 1):
            seq_type = (sorted_type_id[atom] +
                        cumulative_num_types[sorted_elem_id[atom]]
                                            [sorted_spec_id[atom]])
            elem = settings.element_list[sorted_elem_id[atom]]
            f.write(f"{atom:5d} {seq_type:5d} "
                    f"{sorted_direct_xyz[atom][1]:20.14f} "
                    f"{sorted_direct_xyz[atom][2]:20.14f} "
                    f"{sorted_direct_xyz[atom][3]:20.14f} "
                    f"{elem:>2s}\n")

        # Potential sites (identical to atom sites in OLCAO).
        f.write("NUM_POTENTIAL_SITES\n")
        f.write(f"{settings.num_atoms:<5d}\n")
        f.write("NUM_TYPE_X_Y_Z_ELEM\n")
        for atom in range(1, settings.num_atoms + 1):
            seq_type = (sorted_type_id[atom] +
                        cumulative_num_types[sorted_elem_id[atom]]
                                            [sorted_spec_id[atom]])
            elem = settings.element_list[sorted_elem_id[atom]]
            f.write(f"{atom:5d} {seq_type:5d} "
                    f"{sorted_direct_xyz[atom][1]:20.14f} "
                    f"{sorted_direct_xyz[atom][2]:20.14f} "
                    f"{sorted_direct_xyz[atom][3]:20.14f} "
                    f"{elem:>2s}\n")


def _print_olcao_input(settings, sc, file_set,
                       sorted_elem_id, sorted_spec_id, sorted_type_id,
                       sorted_xanes_atoms, cumulative_num_types,
                       total_num_types, num_atoms_of_type):
    """Write ``olcao.dat`` — the main OLCAO control file.

    This is the largest output file produced by makeinput.  It contains
    everything the OLCAO programs need beyond the structure and k-points:

    * Title
    * Atomic orbital basis set for each type (from contracted wave functions)
    * Potential function for each type (SCF starting potential)
    * Exchange-correlation mesh parameters
    * State counts and electron count
    * SCF convergence parameters
    * DOS, bond, SYBD, OPTC, NLOP, SIGE, PACS, field, LOEN parameters
    * END_OF_DATA marker

    Parameters
    ----------
    settings : ScriptSettings
    sc : StructureControl
    file_set : int
    sorted_elem_id, sorted_spec_id, sorted_type_id : lists (1-indexed)
    sorted_xanes_atoms : list
    cumulative_num_types : list of lists (1-indexed)
    total_num_types : int
    num_atoms_of_type : list (1-indexed)

    Returns
    -------
    num_total_core_states : list
        [0, MB, FB, EB] total core states per basis set type.
    num_total_vale_states : list
        [0, MB, FB, EB] total valence states per basis set type.
    num_core_xanes : int
        Number of core orbitals that could be excited (for XANES).
    xanes_control : list
        XANES excitation control strings (1-indexed).
    max_num_atom_alphas : int
        Maximum number of atom Gaussian alphas across all types.
    pot_dim : int
        Total number of potential terms in this file set.
    nuc_charge : float
        Total nuclear charge in the system.
    """
    from math import floor

    olcao_path = os.path.join(INPUTS_DIR, OLCAO_DAT)
    olcao_fh = open(olcao_path, "w")

    # Title.
    olcao_fh.write("TITLE\n")
    title_str = " ".join(sc.system_title[1:])
    olcao_fh.write(f"{title_str}\n")
    olcao_fh.write("END_TITLE\n")

    # Number of atom types.
    olcao_fh.write("NUM_ATOM_TYPES\n")
    olcao_fh.write(f"{total_num_types:<5d}\n")

    # Write each type's basis set and extract state count information.
    (num_total_core_states, num_total_vale_states,
     num_core_xanes, xanes_control, max_num_atom_alphas) = \
        _print_basis_set(
            settings, sc, file_set, olcao_fh,
            sorted_elem_id, sorted_spec_id, sorted_type_id,
            sorted_xanes_atoms, cumulative_num_types,
            total_num_types, num_atoms_of_type)

    # Number of potential types.
    olcao_fh.write("NUM_POTENTIAL_TYPES\n")
    olcao_fh.write(f"{total_num_types:<5d}\n")

    # Write each type's potential function.
    pot_dim, nuc_charge = _print_scf_pot(
        settings, sc, file_set, olcao_fh,
        sorted_elem_id, sorted_spec_id, sorted_type_id,
        sorted_xanes_atoms, cumulative_num_types,
        total_num_types, num_atoms_of_type)

    # Exchange correlation mesh.
    olcao_fh.write("PRINT_XC_MESH\n")
    olcao_fh.write("0\n")
    olcao_fh.write("NUM_ANGULAR_SAMPLE_VECTORS\n")
    olcao_fh.write(f"{settings.num_samp_vectors}\n")
    olcao_fh.write("WTIN_WTOUT\n")
    olcao_fh.write(f"{settings.xc_in_weight} {settings.xc_out_weight}\n")
    olcao_fh.write("RADIAL_SAMPLE-IN_OUT_SPACING\n")
    olcao_fh.write(f"{settings.xc_in_samp} {settings.xc_out_samp}"
                   f" {settings.xc_spacing_samp}\n")

    # Compute the number of electrons and states.
    num_electrons, num_states_used = _get_elec_and_states(
        settings, nuc_charge, num_total_core_states, num_total_vale_states)

    # Shared input data section.
    olcao_fh.write("SHARED_INPUT_DATA\n")

    # Space group.
    olcao_fh.write("SPACE_GROUP\n")
    olcao_fh.write(f"{settings.space_group_num} {settings.space_group_sub_num}\n")

    # Cutoffs.
    olcao_fh.write("BASISFUNCTION_AND_ELECTROSTATIC_CUTOFFS\n")
    olcao_fh.write(f"    1.00000000E-{settings.bf_cutoff}"
                   f"    1.00000000E-{settings.es_cutoff}\n")

    # States.
    olcao_fh.write("NUM_STATES_TO_USE\n")
    olcao_fh.write(f"{num_states_used[1]} {num_states_used[2]}"
                   f" {num_states_used[3]}\n")

    # Electrons.
    olcao_fh.write("NUM_ELECTRONS\n")
    olcao_fh.write(f"{int(num_electrons)}\n")

    # Thermal smearing.
    olcao_fh.write("THERMAL_SMEARING_SIGMA__SIGMA_CUTOFF\n")
    olcao_fh.write(f"{settings.therm_smear_main} 11.5\n")

    # Fermi level search limit.
    olcao_fh.write("FERMI_LEVEL_SEARCH_LIMIT\n")
    olcao_fh.write("13.6\n")

    # Main (SCF) input data.
    olcao_fh.write("MAIN_INPUT_DATA\n")
    olcao_fh.write("LAST_ITERATION\n")
    olcao_fh.write(f"{settings.num_iter_main}\n")
    olcao_fh.write("CONVERGENCE_TEST\n")
    olcao_fh.write(f"{settings.converg_main}\n")
    olcao_fh.write("XC_CODE\n")
    olcao_fh.write(f"{settings.xc_code}\n")
    olcao_fh.write("FEEDBACK_LEVEL\n")
    olcao_fh.write("2\n")
    olcao_fh.write("RELAXATION_FACTOR\n")
    olcao_fh.write("0.2\n")
    olcao_fh.write("EACH_ITER_FLAGS__TDOS\n")
    olcao_fh.write(f"{settings.iter_tdos}\n")
    olcao_fh.write("PLUSUJ_FORM__NUM_ITEMS\n")
    olcao_fh.write("0 0      ! 0,1,2: elements, types, atoms; num items\n")
    olcao_fh.write("ID__U__J\n")
    olcao_fh.write("NUM_SPLIT_TYPES__DEFAULT_SPLIT\n")
    olcao_fh.write(f"0 {settings.spin_main}\n")
    olcao_fh.write("TYPE_ID__SPIN_SPLIT_FACTOR\n")

    # Dipole input data.
    olcao_fh.write("DIPOLE_INPUT_DATA\n")
    olcao_fh.write("DIPOLE_CENTER\n")
    olcao_fh.write("0.0 0.0 0.0    ! FractABC: center about which"
                   " to compute the dipole.\n")

    # DOS input data.
    olcao_fh.write("DOS_INPUT_DATA\n")
    olcao_fh.write(f"{settings.e_delta_dos} {settings.sigma_dos}"
                   f"              ! DOS Delta Energy, DOS sigma broadening"
                   f" (FWHM)\n")
    olcao_fh.write(f"{settings.emin_dos} {settings.emax_dos}"
                   f"                 ! DOS EMIN and EMAX\n")
    olcao_fh.write(f"{settings.detail_code_pdos}"
                   f" ! PDOS ctrl flag: 0 t; 1 a total; 2 a nl; 3 a nlm\n")

    # Bond input data.
    olcao_fh.write("BOND_INPUT_DATA\n")
    olcao_fh.write(f"{settings.max_len_bond}"
                   f"                     ! MAXIMUM BOND LENGTH\n")
    olcao_fh.write(f"{settings.e_delta_bond} {settings.sigma_bond}"
                   f"               ! BOND Delta Energy,"
                   f" BOND sigma broadening\n")
    olcao_fh.write(f"{settings.emin_bond} {settings.emax_bond}"
                   f"                  ! BOND EMIN and EMAX\n")
    olcao_fh.write(f"{settings.output_bondq}"
                   f"        ! Flag for BONDQ output style 0=new a; 1=old t\n")
    olcao_fh.write(f"{settings.bond_3c}"
                   f"                         ! Flag for 3C-BOND\n")
    olcao_fh.write(f"{settings.max_neighbor_bond}"
                   f"                     ! Max # neighbor atoms\n")

    # SYBD input data.
    olcao_fh.write("SYBD_INPUT_DATA\n")
    _process_sybd_path(settings, sc, olcao_fh)

    # MTOP input data.
    olcao_fh.write("MTOP_INPUT_DATA\n")
    olcao_fh.write(f"{settings.kp_mesh_pscf[0]} {settings.kp_mesh_pscf[1]}"
                   f" {settings.kp_mesh_pscf[2]} # PSCF\n")
    olcao_fh.write("0.0 0.0 0.0 # Shift\n")

    # PACS (XANES/ELNES) input data.
    olcao_fh.write("PACS_INPUT_DATA\n")
    olcao_fh.write(f"{sorted_xanes_atoms[file_set]}"
                   f"                       ! Excited atom number\n")
    olcao_fh.write(f"{settings.e_delta_pacs} {settings.sigma_pacs}"
                   f"               ! PACS delta Energy,"
                   f" PACS sigma factor\n")
    olcao_fh.write(f"{settings.onset_slack_pacs}"
                   f" {settings.energy_window_pacs}"
                   f"               ! Energy slack before onset,"
                   f" energy window\n")
    olcao_fh.write(f"{num_core_xanes}"
                   f"                    ! Number of possible core"
                   f" orbitals to excite\n")
    for orbital in range(1, num_core_xanes + 1):
        olcao_fh.write(f"{xanes_control[orbital]}"
                       f" -1 ! QN_n QN_l Init1 Init2 TEDiff\n")

    # Optical properties input data.
    olcao_fh.write("OPTC_INPUT_DATA\n")
    olcao_fh.write(f"{settings.e_cutoff_optc}"
                   f"                    ! OPTC energy cutoff\n")
    olcao_fh.write(f"{settings.e_trans_optc}"
                   f"                     ! OPTC energy trans\n")
    olcao_fh.write(f"{settings.e_delta_optc}"
                   f"                     ! OPTC delta energy\n")
    olcao_fh.write(f"{settings.sigma_optc}"
                   f"                      ! OPTC broadening\n")
    olcao_fh.write(f"{settings.detail_code_poptc}"
                   f"                ! POPTC 0N;1t;2a;3enl;4enlm\n")

    # Nonlinear optical properties.
    olcao_fh.write("NLOP_INPUT_DATA\n")
    olcao_fh.write(f"{settings.e_cutoff_nlop}"
                   f"                    ! NLOP energy cutoff\n")
    olcao_fh.write(f"{settings.e_trans_nlop}"
                   f"                     ! NLOP energy trans\n")
    olcao_fh.write(f"{settings.e_delta_nlop}"
                   f"                     ! NLOP delta energy\n")
    olcao_fh.write(f"{settings.sigma_nlop}"
                   f"                      ! NLOP broadening\n")

    # Sigma(E) input data.
    olcao_fh.write("SIGE_INPUT_DATA\n")
    olcao_fh.write(f"{settings.e_cutoff_sige}"
                   f"                    ! SIGE energy cutoff\n")
    olcao_fh.write(f"{settings.e_trans_sige}"
                   f"                     ! SIGE energy trans\n")
    olcao_fh.write(f"{settings.e_delta_sige}"
                   f"                     ! SIGE delta energy\n")
    olcao_fh.write(f"{settings.sigma_sige}"
                   f"                      ! SIGE broadening\n")

    # Field input data.
    olcao_fh.write("FIELD_INPUT_DATA\n")
    olcao_fh.write("10 10 10"
                   "                      ! a,b,c # of mesh points\n")
    olcao_fh.write("-100000.0 100000.0"
                   "            ! min,max range of energy\n")
    olcao_fh.write("1 1 1 1"
                   "                       ! 0/1 for psi, psi^2, rho, pot\n")
    olcao_fh.write("1"
                   "                             ! 1=Produce XDMF+HDF5 output\n")
    olcao_fh.write("1"
                   "                             ! 1=Produce 1D profile output\n")
    olcao_fh.write("1"
                   "                             ! 1=Produce dipole output\n")

    # Local environment (LOEN) input data.
    olcao_fh.write("LOEN_INPUT_DATA\n")
    olcao_fh.write("1"
                   "                          ! Method: 1 Bispec. Comp.\n")
    olcao_fh.write("4 4"
                   "                        ! Bispec-Comp.: 2j1 2j2\n")
    olcao_fh.write("20 5.0 0.85"
                   "                ! max_neigh cutoff angleSqueeze\n")

    # End of data marker.
    olcao_fh.write("END_OF_DATA\n")

    olcao_fh.close()

    return (num_total_core_states, num_total_vale_states,
            num_core_xanes, xanes_control,
            max_num_atom_alphas, pot_dim, nuc_charge)


def _print_basis_set(settings, sc, file_set, olcao_fh,
                     sorted_elem_id, sorted_spec_id, sorted_type_id,
                     sorted_xanes_atoms, cumulative_num_types,
                     total_num_types, num_atoms_of_type):
    """Write the basis set section of olcao.dat and extract state counts.

    For each (element, species, type) triplet, this function:
    1. Determines the correct basis set file (normal or nocore for XANES).
    2. Writes the ``ATOM_TYPE_ID__SEQUENTIAL_NUMBER`` header and the type
       label.
    3. Appends the entire contracted wave-function file.
    4. Parses the file to extract:
       - The number of Gaussian terms per orbital angular momentum channel.
       - The number of core and valence orbitals.
       - The number of core and valence states for each basis set size
         (minimal, full, extended).
    5. For XANES target atoms, determines which core orbitals could be
       excited.

    Parameters
    ----------
    (see _print_olcao_input for parameter descriptions)

    Returns
    -------
    num_total_core_states, num_total_vale_states : lists [0, MB, FB, EB]
    num_core_xanes : int
    xanes_control : list (1-indexed)
    max_num_atom_alphas : int
    """
    from math import floor

    # If this is not a XANES file set, ensure num_core_xanes starts at 0.
    if file_set == 0:
        num_core_xanes = 0
    else:
        num_core_xanes = 0  # will be updated if needed

    # Initialise total state counts for MB, FB, EB.
    num_total_core_states = [0, 0, 0, 0]
    num_total_vale_states = [0, 0, 0, 0]
    max_num_atom_alphas = settings.max_num_atom_alphas
    xanes_control = [None]  # 1-indexed, filled for XANES

    for element in range(1, settings.num_elements + 1):
        for species in range(1, settings.num_species[element] + 1):
            for type_id in range(1, settings.num_types[file_set][element][species] + 1):
                seq_type = type_id + cumulative_num_types[element][species]

                # Determine if this type is the XANES target atom's type.
                is_xanes = 0
                if (settings.xanes == 1 and file_set > 0 and
                        sorted_xanes_atoms[file_set] > 0):
                    xa = sorted_xanes_atoms[file_set]
                    if (sorted_elem_id[xa] == element and
                            sorted_spec_id[xa] == species and
                            sorted_type_id[xa] == type_id):
                        is_xanes = 1

                # Select the appropriate basis file and label.
                if is_xanes == 1:
                    curr_basis = settings.basis_files_nc[element][species]
                    elem_name = settings.element_list[element]
                    type_label = f"C{elem_name}{species}_{type_id}"
                    incl_core_basis = settings.basis_files[element][species]
                else:
                    curr_basis = settings.basis_files[element][species]
                    elem_name = settings.element_list[element]
                    type_label = f"{elem_name}{species}_{type_id}"

                # Write the header.
                olcao_fh.write("ATOM_TYPE_ID__SEQUENTIAL_NUMBER\n")
                olcao_fh.write(f"{element:<5d} {species:<5d} {type_id:<5d}"
                               f"     {seq_type}\n")
                olcao_fh.write("ATOM_TYPE_LABEL\n")
                olcao_fh.write(f"{type_label}\n")

                # Append the basis file contents.
                basis_path = os.path.join(INPUT_TEMP, curr_basis)
                with open(basis_path, "r") as bf:
                    basis_lines = bf.readlines()
                olcao_fh.writelines(basis_lines)

                # Parse the basis file for state information.
                # Line 0: header line (skip)
                # Line 1: number of terms per orbital angular momentum
                line_idx = 1
                num_terms = [int(x) for x in basis_lines[line_idx].split()]
                num_term_lines = []
                for nt in num_terms:
                    if nt % 4 == 0:
                        num_term_lines.append(nt // 4)
                    else:
                        num_term_lines.append(nt // 4 + 1)

                if num_terms[0] > max_num_atom_alphas:
                    max_num_atom_alphas = num_terms[0]

                # Skip: alphas header (1 line) + alpha values
                line_idx += 1  # now on alphas header
                line_idx += 1  # now on first alpha line
                line_idx += num_term_lines[0]  # past alpha lines

                # Core orbitals header + count
                line_idx += 1  # "NUM_CORE_ORBITALS" or similar header
                core_vals = [int(x) for x in basis_lines[line_idx].split()]
                num_core_orbitals = core_vals[0]  # MB count (same for all)

                # Get core states.
                if num_core_orbitals > 0:
                    num_core_states = _get_orbital_states(
                        basis_lines, line_idx + 1, num_core_orbitals,
                        num_term_lines)
                else:
                    num_core_states = [0, 0, 0, 0]

                # Skip past core orbital data to get to valence.
                # We need to find the valence count line.  After the core
                # orbitals header line, there's an NL_RADIAL_FUNCTIONS header,
                # then the orbital data for each core orbital.
                skip = _count_orbital_lines(
                    basis_lines, line_idx + 1, num_core_orbitals,
                    num_term_lines)
                line_idx += 1 + skip  # past NL_RADIAL_FUNCTIONS header + data

                # Valence orbitals header + count
                line_idx += 1  # header line
                vale_vals = [int(x) for x in basis_lines[line_idx].split()]
                num_vale_orbitals = vale_vals[2] if len(vale_vals) > 2 else vale_vals[0]  # EB

                if num_vale_orbitals > 0:
                    num_vale_states = _get_orbital_states(
                        basis_lines, line_idx + 1, num_vale_orbitals,
                        num_term_lines)
                else:
                    num_vale_states = [0, 0, 0, 0]

                # Accumulate total states.
                for basis_set in range(1, 4):
                    num_total_core_states[basis_set] += (
                        num_core_states[basis_set] *
                        num_atoms_of_type[seq_type])
                    num_total_vale_states[basis_set] += (
                        num_vale_states[basis_set] *
                        num_atoms_of_type[seq_type])

                # For XANES target atoms, determine excitable core orbitals.
                if is_xanes == 1:
                    incl_path = os.path.join(INPUT_TEMP, incl_core_basis)
                    with open(incl_path, "r") as ibf:
                        incl_lines = ibf.readlines()

                    # Parse the included-core file to find core orbital specs.
                    il = 1  # terms line
                    i_terms = [int(x) for x in incl_lines[il].split()]
                    i_term_lines = []
                    for it in i_terms:
                        if it % 4 == 0:
                            i_term_lines.append(it // 4)
                        else:
                            i_term_lines.append(it // 4 + 1)

                    # Skip to core orbital count.
                    il += 1 + 1 + i_term_lines[0]  # header + alphas header + alphas
                    il += 1  # core orbitals header
                    i_core_vals = [int(x) for x in incl_lines[il].split()]
                    num_core_xanes = i_core_vals[0]

                    # Read NL_RADIAL_FUNCTIONS header.
                    il += 1

                    xanes_control = [None]
                    for orbital in range(1, num_core_xanes + 1):
                        # Number of components for this orbital.
                        il += 1
                        comp_vals = [int(x) for x in incl_lines[il].split()]
                        num_components = comp_vals[0]

                        for component in range(1, num_components + 1):
                            il += 1
                            vals = incl_lines[il].split()
                            qn_n = int(vals[0])
                            qn_l = int(vals[1])

                            # Look up initial XANES states.
                            if (qn_n in INIT_XANES_STATES and
                                    qn_l in INIT_XANES_STATES[qn_n]):
                                init1, init2 = INIT_XANES_STATES[qn_n][qn_l]
                            else:
                                init1, init2 = 0, 0

                            xanes_control.append(
                                f"{qn_n} {qn_l} {init1} {init2}")

                            # Skip coefficient lines.
                            if qn_l < len(i_term_lines):
                                il += i_term_lines[qn_l]

    return (num_total_core_states, num_total_vale_states,
            num_core_xanes, xanes_control, max_num_atom_alphas)


def _get_orbital_states(basis_lines, start_line, num_orbitals, num_term_lines):
    """Count the number of states represented by a set of orbitals.

    This parses the NL_RADIAL_FUNCTIONS section of a contracted wave-function
    file.  Each orbital has a "basis tag" (1=MB, 2=FB, 3=EB) indicating
    which basis set sizes include it.  The number of states contributed by
    each orbital component is found in the orbital specification line
    (field index 3), divided by 2 (because the file stores spin-degenerate
    counts, e.g. s states have numStates=2 representing 1 spatial state).

    Parameters
    ----------
    basis_lines : list of str
        The full file contents.
    start_line : int
        Line index of the NL_RADIAL_FUNCTIONS header.
    num_orbitals : int
    num_term_lines : list of int
        Number of coefficient lines per angular momentum channel.

    Returns
    -------
    num_states : list
        [0, MB_states, FB_states, EB_states].
    """
    num_states = [0, 0, 0, 0]
    il = start_line  # NL_RADIAL_FUNCTIONS header
    il += 1  # skip header

    for orbital in range(num_orbitals):
        vals = basis_lines[il].split()
        num_components = int(vals[0])
        basis_tag = int(vals[1])
        il += 1

        for component in range(num_components):
            vals = basis_lines[il].split()
            qn_l = int(vals[1])
            n_states = int(vals[3])
            il += 1

            for basis_set in range(1, 4):
                if basis_tag <= basis_set:
                    num_states[basis_set] += n_states // 2

            # Skip coefficient lines.
            if qn_l < len(num_term_lines):
                il += num_term_lines[qn_l]

    return num_states


def _count_orbital_lines(basis_lines, start_line, num_orbitals,
                         num_term_lines):
    """Count the number of lines occupied by orbital data.

    Used to skip past the core orbital section to reach the valence section.

    Parameters
    ----------
    basis_lines : list of str
    start_line : int
        Line index of the NL_RADIAL_FUNCTIONS header.
    num_orbitals : int
    num_term_lines : list of int

    Returns
    -------
    int
        Total number of lines (including the header).
    """
    il = start_line
    il += 1  # NL_RADIAL_FUNCTIONS header

    for orbital in range(num_orbitals):
        vals = basis_lines[il].split()
        num_components = int(vals[0])
        il += 1

        for component in range(num_components):
            vals = basis_lines[il].split()
            qn_l = int(vals[1])
            il += 1
            if qn_l < len(num_term_lines):
                il += num_term_lines[qn_l]

    return il - start_line


def _print_scf_pot(settings, sc, file_set, olcao_fh,
                   sorted_elem_id, sorted_spec_id, sorted_type_id,
                   sorted_xanes_atoms, cumulative_num_types,
                   total_num_types, num_atoms_of_type):
    """Write the SCF potential section of olcao.dat and the scfV.dat file.

    For each (element, species, type) triplet:
    1. Write the ``POTENTIAL_TYPE_ID__SEQUENTIAL_NUMBER`` header and label
       into ``olcao.dat``.
    2. Read the potential function definition from the database file and
       write it into ``olcao.dat``, applying any potential modifications
       if ``-modpot`` was specified on the command line.
    3. Read the potential coefficients and accumulate them for ``scfV.dat``.
    4. Track the total nuclear charge and potential dimension.

    The ``scfV.dat`` file contains the starting SCF potential in Gaussian
    analytical form.  It has spin-up and spin-down sections (both identical
    to start) and a stub for +U/+J terms.

    Parameters
    ----------
    (see _print_olcao_input)

    Returns
    -------
    pot_dim : int
        Total number of potential Gaussian terms.
    nuc_charge : float
        Total nuclear charge summed over all atoms.
    """
    scf_path = os.path.join(INPUTS_DIR, POTENTIAL)
    scf_fh = open(scf_path, "w")

    scf_fh.write(f"NUM_TYPES {total_num_types}\n")

    nuc_charge = 0.0
    pot_dim = 0
    scf_coeffs = []

    for element in range(1, settings.num_elements + 1):
        for species in range(1, settings.num_species[element] + 1):
            for type_id in range(1, settings.num_types[file_set][element][species] + 1):
                seq_type = type_id + cumulative_num_types[element][species]
                elem_name = settings.element_list[element]

                # Determine XANES label.
                is_xanes = 0
                if (settings.xanes == 1 and file_set > 0 and
                        sorted_xanes_atoms[file_set] > 0):
                    xa = sorted_xanes_atoms[file_set]
                    if (sorted_elem_id[xa] == element and
                            sorted_spec_id[xa] == species and
                            sorted_type_id[xa] == type_id):
                        is_xanes = 1

                if is_xanes == 1:
                    type_label = f"C{elem_name}{species}_{type_id}"
                else:
                    type_label = f"{elem_name}{species}_{type_id}"

                olcao_fh.write("POTENTIAL_TYPE_ID__SEQUENTIAL_NUMBER\n")
                olcao_fh.write(f"{element:<5d} {species:<5d} {type_id:<5d}"
                               f"     {seq_type}\n")
                olcao_fh.write("POTENTIAL_TYPE_LABEL\n")
                olcao_fh.write(f"{type_label}\n")

                # Read the potential function definition.
                pot_file = settings.pot_files[element][species]
                pot_path = os.path.join(INPUT_TEMP, pot_file)
                with open(pot_path, "r") as pf:
                    pot_lines = pf.readlines()

                # Line 0: species header
                olcao_fh.write(pot_lines[0])

                # Line 1: nuclear charge line
                olcao_fh.write(pot_lines[1])
                charge_vals = pot_lines[1].split()
                nuc_charge += float(charge_vals[0]) * num_atoms_of_type[seq_type]

                # Lines 2-4: covalent radius, value, num alphas
                olcao_fh.write(pot_lines[2])
                olcao_fh.write(pot_lines[3])
                olcao_fh.write(pot_lines[4])

                # Lines 5-6-7: num terms, ALPHAS label, alpha range
                # Check for potential modification.
                if (settings.mod_element_name.lower() ==
                        elem_name.lower() and settings.mod_pot == 1):
                    olcao_fh.write(f"{settings.num_mod_terms}\n")
                    pot_dim += settings.num_mod_terms
                    olcao_fh.write("ALPHAS\n")
                    olcao_fh.write(f"{settings.min_mod_term:14.8e}"
                                   f"  {settings.max_mod_term:14.8e}\n")
                else:
                    olcao_fh.write(pot_lines[5])
                    nt_vals = pot_lines[5].split()
                    pot_dim += int(nt_vals[0])
                    olcao_fh.write(pot_lines[6])
                    olcao_fh.write(pot_lines[7])

                # Read the coefficient file.
                coeff_file = settings.coeff_files[element][species]
                coeff_path = os.path.join(INPUT_TEMP, coeff_file)
                with open(coeff_path, "r") as cf:
                    temp_coeffs = cf.readlines()

                # Apply potential modification to coefficients if needed.
                if (settings.mod_element_name.lower() ==
                        elem_name.lower() and settings.mod_pot == 1):
                    n_existing = len(temp_coeffs) - 1  # first line is count
                    if settings.num_mod_terms < n_existing:
                        temp_coeffs = temp_coeffs[:settings.num_mod_terms + 1]
                        temp_coeffs[0] = f" {settings.num_mod_terms}\n"
                    elif settings.num_mod_terms > n_existing:
                        for _ in range(n_existing, settings.num_mod_terms):
                            temp_coeffs.append(
                                "    0.00000000        0.00000000"
                                " 0.0 0.0 0.0\n")
                        temp_coeffs[0] = f" {settings.num_mod_terms}\n"

                scf_coeffs.extend(temp_coeffs)

    # Write spin-up, spin-down, and +U/+J stub.
    scf_fh.write(" TOTAL__OR__SPIN_UP\n")
    scf_fh.writelines(scf_coeffs)
    scf_fh.write(" SPIN_DN\n")
    scf_fh.writelines(scf_coeffs)
    scf_fh.write("NUM__PLUSUJ__TERMS\n")
    scf_fh.write("0\n")
    scf_fh.write(" TOTAL__OR__SPIN_UP\n")
    scf_fh.write(" SPIN_DN\n")

    scf_fh.close()

    return pot_dim, nuc_charge


def _get_elec_and_states(settings, nuc_charge, num_total_core_states,
                         num_total_vale_states):
    """Compute the number of electrons and states to use.

    After orthogonalization, the number of electrons in the calculation is
    the total nuclear charge minus the core-state charge.  The number of
    states for each basis size is the lesser of (num_electrons × state_factor)
    and the total number of valence states.

    Parameters
    ----------
    settings : ScriptSettings
    nuc_charge : float
    num_total_core_states : list [0, MB, FB, EB]
    num_total_vale_states : list [0, MB, FB, EB]

    Returns
    -------
    num_electrons : int
    num_states_used : list [0, MB, FB, EB]
    """
    num_electrons = int(nuc_charge - num_total_core_states[1] * 2)

    num_states_used = [0, 0, 0, 0]
    for basis in range(1, 4):
        if num_electrons * settings.state_factor < num_total_vale_states[basis]:
            if num_electrons % 2 == 0:
                num_states_used[basis] = int(num_electrons * settings.state_factor)
            else:
                num_states_used[basis] = int(num_electrons * settings.state_factor) + 1
        else:
            num_states_used[basis] = num_total_vale_states[basis]

    return num_electrons, num_states_used


def _process_sybd_path(settings, sc, olcao_fh):
    """Process the SYBD high-symmetry path file and write it to olcao.dat.

    The SYBD (Symmetric Band Structure) path file lives in the sybd database
    and defines the high-symmetry k-point path through the Brillouin zone.
    The file uses symbolic variables (lattice constants a, b, c; angles
    alpha, beta, gamma; and derived quantities like eta, zeta) that must be
    substituted with numerical values for the specific crystal.

    Processing steps:

    1. Read the number of lattice-specific variables and their defining
       equations.
    2. Substitute angle names (alpha, beta, gamma) with their radian values.
    3. Substitute references to previously defined variables.
    4. Substitute lattice magnitudes (written as `` a ``, `` b ``, `` c ``
       with surrounding spaces to avoid collisions with substrings of
       function names like cos, csc, etc.).
    5. Evaluate each equation to get a numerical value.
    6. Read the path specification: number of paths, k-points per path,
       k-points per segment, and the high-symmetry k-point coordinates.
    7. For each high-symmetry k-point, substitute variables and evaluate
       the coordinate expressions.
    8. Write everything to olcao.dat.
    9. If ``-printbz`` is active, also append to a BZ.1 file.

    Parameters
    ----------
    settings : ScriptSettings
    sc : StructureControl
    olcao_fh : file handle
        The open olcao.dat file.
    """
    from structure_control import PI

    sybd_path = os.path.join(settings.sybd_db, settings.sybd_path)
    with open(sybd_path, "r") as sf:
        sybd_lines = sf.readlines()

    # Strip and split helper.
    def prep(line):
        return line.strip().split()

    line_idx = 0

    # Read the number of lattice-specific variables.
    vals = prep(sybd_lines[line_idx])
    num_latt_vars = int(vals[0])
    line_idx += 1

    # Angle names for substitution (1-indexed: alpha=1, beta=2, gamma=3).
    angle_name = ["", "alpha", "beta", "gamma"]

    variable_name = [None]  # 1-indexed
    variable_eqn = [None]   # 1-indexed (as strings initially, then floats)

    for var in range(1, num_latt_vars + 1):
        parts = sybd_lines[line_idx].strip().split(" = ", 1)
        variable_name.append(parts[0].strip())
        eqn = parts[1].strip() if len(parts) > 1 else parts[0].strip()
        line_idx += 1

        # Substitute angles (in radians).
        for angle in range(1, 4):
            while angle_name[angle] in eqn:
                eqn = eqn.replace(angle_name[angle],
                                  str(sc.full_cell_angle[angle]))

        # Substitute previously defined variables.
        for prev in range(1, var):
            while variable_name[prev] in eqn:
                eqn = eqn.replace(variable_name[prev],
                                  str(variable_eqn[prev]))

        # Substitute lattice magnitudes.  They appear as " a ", " b ", " c "
        # (with surrounding spaces) to avoid colliding with function names
        # like cos, csc, sec, etc.
        for axis in range(1, 4):
            letter = " " + chr(ord("a") + axis - 1) + " "
            while letter in eqn:
                eqn = eqn.replace(letter, str(sc.full_cell_mag[axis]))

        # Evaluate the equation to a number.
        import math
        # Make math functions available for eval.
        safe_ns = {"__builtins__": {}, "cos": math.cos, "sin": math.sin,
                   "tan": math.tan, "sqrt": math.sqrt, "acos": math.acos,
                   "asin": math.asin, "atan": math.atan, "pi": math.pi,
                   "cot": lambda x: 1.0 / math.tan(x),
                   "sec": lambda x: 1.0 / math.cos(x),
                   "csc": lambda x: 1.0 / math.sin(x)}
        variable_eqn.append(eval(eqn, safe_ns))

    # Read the path specification line.
    vals = prep(sybd_lines[line_idx])
    line_idx += 1
    olcao_fh.write(" ".join(vals) + "\n")
    num_paths = int(vals[0])
    num_path_kp = int(vals[1])

    # If printing BZ, open the file for appending.
    bz_fh = None
    if settings.print_bz > 0:
        bz_fh = open("BZ.1", "a")
        bz_fh.write(f"\nData.num_total_BZ_path_KP = {num_path_kp}\n")
        bz_fh.write(f"\nData.num_high_symmetry_BZ_paths = {num_paths}\n")
        bz_fh.write("\nData.num_high_symmetry_BZ_points_per_path = (")

    # Read the number of high-symmetry k-points per path.
    vals = prep(sybd_lines[line_idx])
    line_idx += 1
    olcao_fh.write(" ".join(vals) + "\n")
    num_high_symm_kp = 0
    for i in range(num_paths):
        num_high_symm_kp += int(vals[i])
        if bz_fh is not None:
            if i < num_paths - 1:
                bz_fh.write(f"{vals[i]},")
            else:
                bz_fh.write(f"{vals[i]})\n")

    if bz_fh is not None:
        bz_fh.write("\nData.high_symmetry_BZ_kpoints = [\n")

    # Read and process each high-symmetry k-point.
    symm_kp_char = ""
    for kp in range(1, num_high_symm_kp + 1):
        vals = prep(sybd_lines[line_idx])
        line_idx += 1

        # Substitute variables into the first three coordinate fields.
        for axis in range(3):
            expr = vals[axis]
            for var in range(1, num_latt_vars + 1):
                while variable_name[var] in expr:
                    expr = expr.replace(variable_name[var],
                                        str(variable_eqn[var]))
            # Remove double negatives.
            expr = expr.replace("--", "+")
            vals[axis] = str(eval(expr))

        # Write to olcao.dat.
        olcao_fh.write(f"{float(vals[0]):12.8f}{float(vals[1]):12.8f}"
                       f"{float(vals[2]):12.8f} {vals[3]} {vals[4]:<5s}\n")

        # Collect the first character of the label for the summary string.
        if len(vals) > 4:
            symm_kp_char += vals[4][0]

        # If printing BZ, add these k-points.
        if bz_fh is not None:
            if kp < num_high_symm_kp:
                bz_fh.write(f"[{float(vals[0]):12.8f} ,"
                            f" {float(vals[1]):12.8f} ,"
                            f" {float(vals[2]):12.8f}],\n")
            else:
                bz_fh.write(f"[{float(vals[0]):12.8f} ,"
                            f" {float(vals[1]):12.8f} ,"
                            f" {float(vals[2]):12.8f}]]\n")

    if bz_fh is not None:
        bz_fh.close()

    # Write the string of high-symmetry k-point label characters.
    olcao_fh.write(f"{symm_kp_char}\n")


def _get_kp():
    """Copy k-point files from .inputTemp into the inputs directory."""
    for kp_file in [KP_SCF, KP_PSCF]:
        src = os.path.join(INPUT_TEMP, kp_file)
        dst = os.path.join(INPUTS_DIR, kp_file)
        if os.path.exists(src):
            shutil.copy2(src, dst)


def _make_olcao_mi(settings, sc, sorted_elem_id, sorted_spec_id,
                   sorted_fract_abc, sorted_direct_xyz, output_file):
    """Write an olcao.mi skeleton file with explicit species assignments.

    This function creates a new skeleton file (either fractional or Cartesian)
    where each atom's species assignment is explicitly recorded.  The file
    can be used to re-run makeinput with the same species assignments without
    relying on the grouping methods.

    Two versions are created:
    * ``olcao.fract-mi`` — atomic positions in fractional coordinates
    * ``olcao.cart-mi``  — atomic positions in Cartesian coordinates (Å)

    Parameters
    ----------
    settings : ScriptSettings
    sc : StructureControl
    sorted_elem_id, sorted_spec_id : lists (1-indexed)
    sorted_fract_abc, sorted_direct_xyz : lists of lists (1-indexed)
    output_file : str
        Either OLCAO_FMI or OLCAO_CMI.
    """
    from structure_control import BOHR_RAD, PI

    path = os.path.join(INPUTS_DIR, output_file)
    with open(path, "w") as f:
        f.write("title\n")
        f.write("OLCAO file created by makeinput\n")
        title_str = " ".join(sc.system_title[1:])
        f.write(f"{title_str}\n")
        f.write("end\n")
        f.write("cell\n")

        # Print cell parameters in a, b, c, alpha, beta, gamma format.
        # Convert magnitudes back to Angstroms (they were converted to Bohr
        # by _convert_a_to_au).
        f.write(f"{sc.mag[1] * BOHR_RAD:13.8f} "
                f"{sc.mag[2] * BOHR_RAD:13.8f} "
                f"{sc.mag[3] * BOHR_RAD:13.8f} "
                f"{sc.angle[1] * 180.0 / PI:13.8f} "
                f"{sc.angle[2] * 180.0 / PI:13.8f} "
                f"{sc.angle[3] * 180.0 / PI:13.8f}\n")

        # Atom positions header.
        if output_file == OLCAO_FMI:
            f.write(f"fract {settings.num_atoms} 1\n")
        else:
            f.write(f"cart {settings.num_atoms} 1\n")

        # Atom positions.
        for atom in range(1, settings.num_atoms + 1):
            elem = settings.element_list[sorted_elem_id[atom]]
            spec = sorted_spec_id[atom]
            if output_file == OLCAO_FMI:
                f.write(f"{elem:>2s}{spec:<4d} "
                        f"{sorted_fract_abc[atom][1]:12.6f} "
                        f"{sorted_fract_abc[atom][2]:12.6f} "
                        f"{sorted_fract_abc[atom][3]:12.6f}\n")
            else:
                f.write(f"{elem:>2s}{spec:<4d} "
                        f"{sorted_direct_xyz[atom][1] * BOHR_RAD:12.6f} "
                        f"{sorted_direct_xyz[atom][2] * BOHR_RAD:12.6f} "
                        f"{sorted_direct_xyz[atom][3] * BOHR_RAD:12.6f}\n")

        # Footer.
        f.write("space 1_a\n")
        f.write("supercell 1 1 1\n")
        f.write("full\n")


def _make_pdb(sc, file_set, current_name, sorted_elem_id, sorted_spec_id,
              sorted_fract_abc, sorted_direct_xyz):
    """Produce a PDB file for structure visualisation.

    This delegates to the StructureControl ``print_pdb`` method, using
    the sorted atom data for this file set.

    Parameters
    ----------
    sc : StructureControl
    file_set : int
    current_name : str
    sorted_elem_id, sorted_spec_id : lists (1-indexed)
    sorted_fract_abc, sorted_direct_xyz : lists of lists (1-indexed)
    """
    sc.print_pdb(filename=current_name + PDB_EXT)


def _make_cif(sc, settings, file_set, current_name, sorted_elem_id,
              sorted_spec_id, sorted_type_id, sorted_fract_abc):
    """Produce a CIF file for structure visualisation.

    Delegates to the StructureControl ``print_cif`` method.

    Parameters
    ----------
    sc : StructureControl
    settings : ScriptSettings
    file_set : int
    current_name : str
    sorted_elem_id, sorted_spec_id, sorted_type_id : lists (1-indexed)
    sorted_fract_abc : list of lists (1-indexed)
    """
    sc.print_cif(filename=current_name + CIF_EXT)


def _compute_and_print_mem(settings, sc, file_set,
                           num_total_core_states, num_total_vale_states,
                           total_num_types, num_atoms_of_type,
                           max_num_atom_alphas, pot_dim, nuc_charge):
    """Compute memory estimates for each OLCAO phase and write them.

    For each of the three basis set sizes (minimal, full, extended), the
    memory usage is estimated for every OLCAO program phase:

    * **setup**: The dominant allocation is in atomPotOverlap — matrices
      whose sizes depend on core/valence state counts and k-point counts.
      The reciprocal-space electrostatic contribution may dominate for large
      cells.
    * **main**: Peak is either during Hamiltonian evaluation (overlap and
      Hamiltonian matrices) or during SCF potential construction (potDim²).
    * **intg**: Fixed estimate based on the hard-coded dataChunkSize.
    * **band**: Overlap + Hamiltonian + orthogonalisation matrices.
    * **PDOS**: Wave-function matrices + output spectra storage.
    * **bond**: Wave-function + bond-order matrices.
    * **OPTC**: Transition pairs + momentum matrices + wave functions.

    The estimates count floating-point numbers, then convert to bytes
    (×OLCAO_DOUBLE) and express in megabytes or gigabytes depending on
    magnitude.

    Parameters
    ----------
    settings : ScriptSettings
    sc : StructureControl
    file_set : int
    num_total_core_states, num_total_vale_states : lists [0, MB, FB, EB]
    total_num_types : int
    num_atoms_of_type : list (1-indexed)
    max_num_atom_alphas : int
    pot_dim : int
    nuc_charge : float
    """
    from math import sqrt, ceil, log

    from structure_control import BOHR_RAD

    # Compute electrons and states (needed for memory estimates).
    num_electrons, num_states_used = _get_elec_and_states(
        settings, nuc_charge, num_total_core_states, num_total_vale_states)

    # Storage for memory estimates: [basis_type] for each phase.
    mem_setup = [0.0, 0.0, 0.0, 0.0]
    mem_main = [0.0, 0.0, 0.0, 0.0]
    mem_intg = [0.0, 0.0, 0.0, 0.0]
    mem_band = [0.0, 0.0, 0.0, 0.0]
    mem_pdos = [0.0, 0.0, 0.0, 0.0]
    mem_bond = [0.0, 0.0, 0.0, 0.0]
    mem_optc = [0.0, 0.0, 0.0, 0.0]

    for bt in range(1, 4):
        nc = num_total_core_states[bt]
        nv = num_total_vale_states[bt]
        kp1 = settings.kp_num[1]  # SCF kpoints
        kp2 = settings.kp_num[2]  # PSCF kpoints

        # ---- setup ----
        mem_s = 0.0
        # alphaPotDist
        mem_s += (settings.num_elements ** 3 *
                  max_num_atom_alphas ** 2 *
                  settings.max_num_pot_alphas)
        # coreValeOL + valeVale + coreCore
        mem_s += nc * nv * kp1
        mem_s += nv * nv * kp1
        mem_s += nc * nc * kp1
        # Compare coreVale*kp vs (valeCore + packed valeVale)
        alt = nc * nv * kp1
        other = nv * nc + nv * (nv + 1) / 2.0
        if alt < other:
            mem_s += nv * nc + nv * (nv + 1) / 2.0
        else:
            mem_s += nc * nv * kp1
        if settings.kp_note[1] == "(General)":
            mem_s *= 2.0

        # Reciprocal-space electrostatic contribution.
        # Convert reciprocal lattice to atomic units for the estimate.
        recip = [[None, 0.0, 0.0, 0.0] for _ in range(4)]
        recip[0] = None
        for i in range(1, 4):
            for j in range(1, 4):
                recip[j][i] = sc.recip_lattice[j][i] * BOHR_RAD

        neglig_limit = 4.0 * -log(10 ** (settings.es_cutoff * -1)) * settings.min_pot_alpha
        prim_reps = [0, 0, 0, 0]
        for i in range(1, 4):
            s = 0.0
            for j in range(1, 4):
                s += recip[j][i] ** 2
            s = sqrt(s)
            prim_reps[i] = 2 * int(ceil(sqrt(neglig_limit) / s)) + 1

        num_cells_recip = prim_reps[1] * prim_reps[2] * prim_reps[3] / 2.0
        mem_alt = num_cells_recip * (settings.max_num_pot_alphas + 4.0)
        if mem_alt > mem_s:
            mem_s = mem_alt
        mem_setup[bt] = mem_s

        # ---- main ----
        mem_m = 0.0
        mem_m += nv * kp1  # energy eigenvalues
        mem_m += nv * nv    # valeValeOL
        mem_m += nv * nv    # valeVale
        mem_m += nv * (nv + 1) / 2  # packedValeVale
        mem_m += nv * (nv + 1) / 2  # tempPackedValeVale
        if settings.kp_note[1] == "(General)":
            mem_m *= 2.0
        # SCF potential peak
        mem_m_alt = (4 * pot_dim * pot_dim +
                     pot_dim * 28 + (pot_dim + 4) * 4000)
        if mem_m_alt > mem_m:
            mem_m = mem_m_alt
        mem_main[bt] = mem_m

        # ---- intg ----
        mem_intg[bt] = 10000000 * 5

        # ---- band ----
        mem_b = (nv * nv * 2 + nc * nv * 2 + nc * nc)
        if settings.kp_note[2] == "(General)":
            mem_b *= 2.0
        mem_b += 10000000
        mem_band[bt] = mem_b

        # ---- PDOS ----
        mem_p = nv * nv + nv * num_states_used[bt]
        if settings.kp_note[2] == "(Gamma)":
            pass
        else:
            mem_p *= 3
        num_energy_pts = (settings.emax_dos - settings.emin_dos) / settings.e_delta_dos
        mem_p += num_energy_pts * (3 + total_num_types)
        mem_p += settings.max_num_vale_states * settings.num_atoms
        mem_pdos[bt] = mem_p

        # ---- bond ----
        mem_bo = nv * nv + nv * num_states_used[bt]
        if settings.kp_note[2] == "(Gamma)":
            pass
        else:
            mem_bo *= 3
        mem_bo += settings.num_atoms * settings.num_atoms * 2
        num_energy_pts_b = (settings.emax_dos - settings.emin_dos) / settings.e_delta_dos
        if file_set > 0:
            mem_bo += num_energy_pts_b * 31
        mem_bond[bt] = mem_bo

        # ---- OPTC ----
        if file_set > 0:
            est_pairs = 5 * 500
        else:
            est_pairs = 500 * 500
        mem_o = est_pairs * kp2 * 4
        mem_wf = nv * num_states_used[bt]
        if settings.kp_note[2] == "(Gamma)":
            mem_o += mem_wf * 2
        else:
            mem_o += mem_wf * 3
        mem_o += nv * nv * 3
        mem_o += nv * nv * 3 * 2
        mem_half = nv * num_states_used[bt] / 2.0
        if settings.kp_note[2] == "(Gamma)":
            mem_o += mem_half * 3
        else:
            mem_o += mem_half * 3 * 2
        mem_optc[bt] = mem_o

    # ---- Adjust to bytes and choose units ----
    all_phases = [mem_setup, mem_main, mem_intg, mem_band,
                  mem_pdos, mem_bond, mem_optc]

    # Add executable overhead.
    for phase in all_phases:
        for bt in range(1, 4):
            phase[bt] += 10000000

    # Find the maximum and choose units.
    max_est = 0.0
    for phase in all_phases:
        for bt in range(1, 4):
            if phase[bt] > max_est:
                max_est = phase[bt]

    max_est_bytes = max_est * OLCAO_DOUBLE
    if max_est_bytes < 100000000:
        memory_note = "Megabytes"
        memory_factor = 1000000
    else:
        memory_note = "Gigabytes"
        memory_factor = 1000000000

    # Convert to the chosen units.
    for phase in all_phases:
        for bt in range(1, 4):
            phase[bt] *= OLCAO_DOUBLE / memory_factor

    # Write the memory file.
    mem_path = os.path.join(INPUTS_DIR, MEMORY_FILE)
    with open(mem_path, "w") as f:
        f.write(f"Memory Estimates:  ({memory_note})\n")
        f.write("    SCF/PSCF KPoints:\n")
        f.write("     SETUP      MAIN      INTG      BAND"
                "      PDOS      BOND      OPTC\n")
        for bt in range(1, 4):
            f.write(f"   {mem_setup[bt]:7.3f}    {mem_main[bt]:7.3f}"
                    f"   {mem_intg[bt]:7.3f}   {mem_band[bt]:7.3f}"
                    f"   {mem_pdos[bt]:7.3f}   {mem_bond[bt]:7.3f}"
                    f"   {mem_optc[bt]:7.3f}\n")


def _make_sub_file(settings, h, current_name, sub_dir, proj_home,
                   num_core_xanes, xanes_control):
    """Write the job submission script (PBS, LSF, Slurm, or bash).

    The submission file is placed in ``sub_dir/`` (which is either ``.``
    for non-XANES or the per-target directory for XANES).  The file
    contains:
    * A system-specific header (scheduler directives or bash shebang).
    * Source of the OLCAO environment.
    * The appropriate OLCAO commands for this file set.

    For XANES file sets (h > 0), the commands iterate over each excitable
    core orbital and run the PACS calculation with the appropriate edge
    label (1s, 2p, 3d, etc.).

    Parameters
    ----------
    settings : ScriptSettings
    h : int
        File set index (0 = base, >0 = XANES target).
    current_name : str
        Name of this file set / directory.
    sub_dir : str
        Directory where the submission file should be written.
    proj_home : str
        Basename of the project directory.
    num_core_xanes : int
    xanes_control : list (1-indexed)
    """
    olcao_rc = os.getenv("OLCAO_RC", "")

    # Define the project name for the queue (truncated to 15 chars).
    if h == 0:
        proj_name = proj_home
    else:
        proj_name = f"{proj_home}-{current_name.lower()}"
    if len(proj_name) > 15:
        proj_name = proj_name[:15]

    proj_dir = os.getcwd()

    # Determine which submission file and header to use.
    if settings.queue_type == 1:
        sub_path = os.path.join(sub_dir, PBS_SUB)
        header = [
            f"#PBS -N {proj_name}\n",
            "#PBS -l cput=100:0:0,ncpus=1\n",
            "#\n",
        ]
    elif settings.queue_type == 2:
        sub_path = os.path.join(sub_dir, LSF_SUB)
        header = [
            f"#BSUB -J {proj_name}\n",
            f"#BSUB -oo {proj_name}.o%J\n",
            f"#BSUB -eo {proj_name}.e%J\n",
            "#BSUB -c 100:0\n",
            "#\n",
        ]
    elif settings.queue_type == 3:
        sub_path = os.path.join(sub_dir, SLURM_SUB)
        header = [
            "#!/bin/bash\n",
            f"#SBATCH -p {settings.partition}\n",
        ]
        if settings.account != "":
            header.append(f"#SBATCH -A {settings.account}\n")
        header.extend([
            f"#SBATCH -J {proj_name}\n",
            f"#SBATCH -o {proj_name}.o%J\n",
            f"#SBATCH -e {proj_name}.e%J\n",
            f"#SBATCH -n {settings.cpus}\n",
            f"#SBATCH -N {settings.nodes}\n",
            f"#SBATCH -t {settings.time}\n",
            f"#SBATCH --mem={settings.memory}\n",
            "#\n",
        ])
    else:
        sub_path = os.path.join(sub_dir, BASH_SUB)
        bash_loc = shutil.which("bash") or "/bin/bash"
        header = [f"#!{bash_loc}\n"]

    with open(sub_path, "w") as f:
        f.writelines(header)

        if settings.xanes == 0:
            # Non-XANES: simple run commands.
            f.write("# shellcheck source=/dev/null\n")
            f.write(f"source {olcao_rc}/olcaorc\n")
            f.write("export OMP_NUM_THREADS=1\n")
            f.write(f"cd {proj_dir} || exit\n")
            f.write("\"$OLCAO_BIN/uolcao\" -bond\n")
        elif h == 0:
            # XANES base: runs in the project subdirectory.
            f.write("# shellcheck source=/dev/null\n")
            f.write(f"source {olcao_rc}/olcaorc\n")
            basename = "/" + os.path.basename(proj_dir)
            f.write("export OMP_NUM_THREADS=1\n")
            f.write(f"cd {proj_dir}{basename} || exit\n")
            f.write("\"$OLCAO_BIN/uolcao\" -dos\n")
            f.write("\"$OLCAO_BIN/uolcao\" -bond\n")
            f.write("\"$OLCAO_BIN/uolcao\" -sybd\n")
            f.write("\"$OLCAO_BIN/uolcao\" -optc\n")
        else:
            # XANES target: PACS commands for each excitable core orbital.
            f.write("# shellcheck source=/dev/null\n")
            f.write(f"source {olcao_rc}/olcaorc\n")
            f.write("export OMP_NUM_THREADS=1\n")
            f.write("cd $SLURM_SUBMIT_DIR || exit\n")

            l_to_char = {0: "s", 1: "p", 2: "d", 3: "f"}
            for orbital in range(1, num_core_xanes + 1):
                vals = xanes_control[orbital].split()
                qn_n = vals[0]
                qn_l = int(vals[1])
                edge_char = l_to_char.get(qn_l, "?")
                f.write(f"\"$OLCAO_BIN/olcao\" -pacs {qn_n}{edge_char}\n")

    # If bash was used, make the submission file executable.
    if settings.queue_type == 0:
        os.chmod(sub_path, os.stat(sub_path).st_mode | 0o755)


def print_summary(settings, sc):
    """Write a concise summary of the generated OLCAO input.

    The summary file provides a quick-reference snapshot of the key
    parameters and dimensions of the calculation that was just set up.
    It is intended to be human-readable so that a user (or student) can
    verify at a glance that the inputs are reasonable before submitting
    a job.

    Contents of the summary file:

    * **Lattice parameters** — a, b, c in Angstroms and alpha, beta, gamma
      in degrees, followed by the real-space lattice vectors in atomic units
      (Bohr).
    * **Matrix dimensions** — the number of valence states and core states
      for each basis set size (minimal, full, extended), plus the total
      number of potential Gaussian terms.  These determine the sizes of the
      matrices that OLCAO will allocate and thus the memory and CPU cost.
    * **System statistics** — number of atoms, electrons, and states to
      use for each basis.
    * **K-point information** — the number of SCF and post-SCF k-points,
      with the mesh dimensions and a note indicating whether a gamma-only
      or general k-point set is being used.
    * **Execution parameters** — convergence limit, thermal smearing, and
      broadening sigmas for DOS, OPTC, and SIGE.
    * **Type count** — the total number of distinct types for the base
      (non-XANES) file set.
    * **XANES details** (if applicable) — for each target atom, the
      potential dimension, type count, and the mapping between sorted
      (DAT) and original (skeleton) atom numbers.

    Parameters
    ----------
    settings : ScriptSettings
        Must have ``_summary_data`` and ``_proj_home`` attributes set by
        a prior call to ``print_olcao``.
    sc : StructureControl
        The structure control instance (lattice data is read from here).
    """
    from structure_control import PI, BOHR_RAD

    # Retrieve the base file set data (h=0).
    sd0 = settings._summary_data[0]
    num_total_vale_states = sd0['num_total_vale_states']
    num_total_core_states = sd0['num_total_core_states']
    num_electrons = sd0['num_electrons']
    num_states_used = sd0['num_states_used']
    proj_home = settings._proj_home

    num_xanes_atoms = len(settings.xanes_atoms)

    # Determine the spacing for k-point mesh formatting.
    kp_index_space = 1
    for axis in range(3):
        if settings.kp_mesh_scf[axis] >= 10:
            kp_index_space = 2
        if settings.kp_mesh_pscf[axis] >= 10:
            kp_index_space = 2

    with open(SUMMARY_FILE, "w") as f:
        # Lattice parameters in Angstroms and degrees.
        f.write(f"    a = {sc.mag[1] * BOHR_RAD:12.8f} A;"
                f"       b = {sc.mag[2] * BOHR_RAD:12.8f} A;"
                f"        c = {sc.mag[3] * BOHR_RAD:12.8f} A\n")
        f.write(f"alpha = {sc.angle[1] * 180.0 / PI:12.8f} Deg.;"
                f" beta = {sc.angle[2] * 180.0 / PI:12.8f} Deg.;"
                f" gamma = {sc.angle[3] * 180.0 / PI:12.8f} Deg.\n")

        # Lattice vectors in atomic units.
        f.write(f"{sc.real_lattice[1][1]:16.8f}"
                f"{sc.real_lattice[1][2]:16.8f}"
                f"{sc.real_lattice[1][3]:16.8f} ax ay az a.u.\n")
        f.write(f"{sc.real_lattice[2][1]:16.8f}"
                f"{sc.real_lattice[2][2]:16.8f}"
                f"{sc.real_lattice[2][3]:16.8f} bx by bz a.u.\n")
        f.write(f"{sc.real_lattice[3][1]:16.8f}"
                f"{sc.real_lattice[3][2]:16.8f}"
                f"{sc.real_lattice[3][3]:16.8f} cx cy cz a.u.\n")

        # Matrix dimensions.
        f.write(f"Dimensions: Vale        ="
                f" {num_total_vale_states[1]}"
                f" {num_total_vale_states[2]}"
                f" {num_total_vale_states[3]}\n")
        f.write(f"Dimensions: Core        ="
                f" {num_total_core_states[1]}"
                f" {num_total_core_states[2]}"
                f" {num_total_core_states[3]}\n")
        f.write(f"Dimension:  Pot         =  {sd0['pot_dim']}\n")

        # System statistics.
        f.write(f"Number of Atoms         =  {settings.num_atoms}\n")
        f.write(f"Number of Electrons     =  {int(num_electrons)}\n")
        f.write(f"Number of States        ="
                f" {num_states_used[1]}"
                f" {num_states_used[2]}"
                f" {num_states_used[3]}\n")

        # K-point information.  In density mode the exact
        # count and mesh dimensions are unknown at this stage
        # (they are computed by OLCAO at runtime), so we print
        # the density value instead.
        intg_names = {0: "Gaussian", 1: "LAT"}
        scf_intg = intg_names.get(
            settings.kp_intg_code[1],
            str(settings.kp_intg_code[1]))
        pscf_intg = intg_names.get(
            settings.kp_intg_code[2],
            str(settings.kp_intg_code[2]))
        if settings.use_kp_density:
            f.write(
                f"SCF KPoints            =  "
                f"density={settings.kp_density[1]}"
                f"  {settings.kp_note[1]}\n")
            f.write(
                f"PSCF KPoints           =  "
                f"density={settings.kp_density[2]}"
                f"  {settings.kp_note[2]}\n")
        else:
            w = kp_index_space
            f.write(
                f"Number of SCF KPoints  =  "
                f"{settings.kp_num[1]:<4d}"
                f" {settings.kp_note[1]:>9s}  "
                f"[{settings.kp_mesh_scf[0]:>{w}d}"
                f" {settings.kp_mesh_scf[1]:>{w}d}"
                f" {settings.kp_mesh_scf[2]:>{w}d}]\n")
            f.write(
                f"Number of PSCF KPoints =  "
                f"{settings.kp_num[2]:<4d}"
                f" {settings.kp_note[2]:>9s}  "
                f"[{settings.kp_mesh_pscf[0]:>{w}d}"
                f" {settings.kp_mesh_pscf[1]:>{w}d}"
                f" {settings.kp_mesh_pscf[2]:>{w}d}]\n")
        f.write(
            f"SCF KP Integration     =  "
            f"{scf_intg}\n")
        f.write(
            f"PSCF KP Integration    =  "
            f"{pscf_intg}\n")

        # Execution parameters.
        f.write(f"Convergence Limit       =  {settings.converg_main}\n")
        f.write(f"Thermal Sigma           =  {settings.therm_smear_main}\n")
        f.write(f"DOS  Sigma              =  {settings.sigma_dos}\n")
        f.write(f"OPTC Sigma              =  {settings.sigma_optc}\n")
        f.write(f"SIGE Sigma              =  {settings.sigma_sige}\n")

        # Type count for the base file set.
        f.write(f"Number of Types         =  {sd0['total_num_types']}"
                f"   (For the {proj_home} case)\n")

        # XANES details.
        if num_xanes_atoms > 0:
            # The self_min_dist is the minimum image distance for the
            # first XANES target atom.  This tells the user the distance
            # at which the atom "sees" its own periodic image, which is
            # relevant for choosing the XANES radius.
            xa1 = settings.xanes_atoms[0]
            if hasattr(sc, 'self_min_dist') and xa1 < len(sc.self_min_dist):
                f.write(f"All target atoms mirror at"
                        f" {sc.self_min_dist[xa1]} A\n")

        for file_set in range(1, num_xanes_atoms + 1):
            sd = settings._summary_data[file_set]
            xa = settings.xanes_atoms[file_set - 1]
            elem = settings.element_list[settings.atom_element_id[xa]]
            spec = settings.atom_species_id[xa]
            f.write("-------------------------------------------------------\n")
            f.write(f"Dimension:  Pot      = {sd['pot_dim']}\n")
            f.write(f"Number of types      = {sd['total_num_types']}"
                    f" (For the {elem}{spec} case)\n")
            f.write(f"DAT# {sd['sorted_xanes_atoms'][file_set]};"
                    f"  SKELETON# {xa}\n")
            f.write("-------------------------------------------------------\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("\n\nmakeinput.py script executing.\n")

    # Parse command line and load defaults.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(f"\nStarting environment setup at................{ts}.")
    settings = ScriptSettings()

    # Set up the execution environment.
    setup_environment(settings)

    # Create the StructureControl instance that will hold the structure.
    from structure_control import StructureControl
    sc = StructureControl()

    # Read the skeleton, apply space group/supercell, init default species.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(f"\nStarting structure initialization at.........{ts}.")
    initialize_cell(settings, sc)

    # Assign species.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(f"\nStarting to assign species at................{ts}.")
    assign_group(settings, sc, "species")

    # Assign types.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(f"\nStarting to assign general types at..........{ts}.")
    assign_group(settings, sc, "types")

    # XANES types (if requested).
    if settings.xanes == 1:
        ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
        print(f"\nStarting to assign XANES types at............{ts}.")
        assign_xanes_types(settings, sc)

    # EMU configuration (if requested).
    if settings.emu == 1:
        ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
        print(f"\nStarting to print emu configuration at.......{ts}.")
        initialize_emu(settings, sc)

    # Write all output files.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(f"\nStarting to print input files at.............{ts}.")
    print_olcao(settings, sc)

    # Summary.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(f"\nStarting to print summary at.................{ts}.")
    print_summary(settings, sc)

    # Done.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(f"\nmakeinput.py script complete at...............{ts}.\n")


if __name__ == "__main__":
    main()
