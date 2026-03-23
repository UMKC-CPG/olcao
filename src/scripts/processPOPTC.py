#!/usr/bin/env python3

"""processPOPTC — Partial Optical Properties post-processor.

PROGRAM:  processPOPTC.py
PURPOSE:  This program will read the raw Epsilon 2 pOptc data
   file with the purpose of collecting label information and
   data along with creating a control file for makePDOS so
   that the raw data for one partial pair can be formatted
   to run OLCAOkkc. The purpose of this script is to allow
   for the calculation of partial ELF, Eps1, Eps1i and ref.

AUTHOR:  Alysse Weigand (Perl original); Python port by Claude.
LAST MODIFIED:  Mar. 20, 2026

USAGE:
   processPOPTC.py [-s SPIN] [-h | --help]

   The -s option selects the spin channel:
      1 = spin up or total (default)
      2 = spin down

When the OLCAO poptc calculation runs, it will produce a
   "raw" structured output file for all partial epsilon2 sets
   in fort.250. This program will read each partial epsilon2
   from the raw file. One at a time, the program will create
   a new file (fort.50) for each partial epsilon2, then it
   will run OLCAOkkc on that file to produce the partial
   epsilon1, ELF, conductivity, etc. Those data sets are then
   appended into new "raw" files that (when complete) will
   contain all partial data for their respective spectral type.

FILE NUMBERING CONVENTION:
   Input files:
      fort.250 — Raw eps2 data (spin 1 / up / total)
      fort.251 — Raw eps2 data (spin 2 / down)
      fort.209 — KKC control data (scaling factors)

   Intermediate files (created and deleted per sequence):
      fort.700 — Control file for makePDOS
      fort.450 — Temp eps2 (input to OLCAOkkc)
      fort.500 — Temp eps1+eps2+ELF (output from OLCAOkkc)
      fort.510 — Temp eps1 total, x, y, z (from OLCAOkkc)
      fort.520 — Temp ELF total, x, y, z (from OLCAOkkc)
      fort.530 — Temp ref. idx. total, x, y, z (from OLCAOkkc)
      fort.540 — Temp ext. coeff. total, x, y, z (from OLCAOkkc)
      fort.550 — Temp eps1i total, x, y, z (from OLCAOkkc)
      fort.560 — Temp reflectivity total, x, y, z (from OLCAOkkc)
      fort.570 — Temp absorp. coeff. total, x, y, z (from OLCAOkkc)

   Output files (spin 1 / spin 2):
      fort.300/301 — Raw poptc eps1, eps2, ELF
      fort.310/311 — Raw poptc eps1 total, x, y, z
      fort.320/321 — Raw poptc ELF total, x, y, z
      fort.330/331 — Raw poptc ref. idx. total, x, y, z
      fort.340/341 — Raw poptc ext. coeff. total, x, y, z
      fort.350/351 — Raw poptc eps1i total, x, y, z
      fort.360/361 — Raw poptc reflectivity total, x, y, z
      fort.370/371 — Raw poptc absorp. coeff. total, x, y, z
"""

import argparse as ap
import os
import subprocess
import sys
from datetime import datetime


# ============================================================
# Utility: prep_line
# ============================================================

def prep_line(file_handle, splitter=None):
    """Read one line, strip it, and split on *splitter*.

    This is the Python equivalent of the Perl
    StructureControl::prepLine subroutine. It reads a single
    line from *file_handle*, strips leading/trailing
    whitespace, and splits on *splitter*. If *splitter* is
    None the line is split on whitespace (the default
    behaviour of str.split()). A leading empty token produced
    by a leading delimiter is removed.

    Parameters
    ----------
    file_handle : file object
        An open file from which to read one line.
    splitter : str or None
        The delimiter string (or regex-style literal) to
        split on. None means split on whitespace.

    Returns
    -------
    list of str
        The tokens produced by the split, with a leading
        empty string removed if present.
    """

    line = file_handle.readline()
    if not line:
        return []

    line = line.strip()

    if splitter is None:
        values = line.split()
    else:
        values = line.split(splitter)

    # Remove a leading empty token that can result from a
    # leading delimiter (mirrors the Perl shift behaviour).
    if values and values[0] == '':
        values.pop(0)

    # Also strip individual tokens when splitting on a
    # non-whitespace delimiter (e.g. ':').
    if splitter is not None:
        values = [v.strip() for v in values]

    return values


# ============================================================
# ScriptSettings: command-line and rc-file handling
# ============================================================

class ScriptSettings():
    """Holds all user-configurable settings for processPOPTC.

    The instance variables of this object are the user
    settings that control the program. The variable values are
    pulled from a list that is created within a resource
    control file and that are then reconciled with command line
    parameters.
    """

    def __init__(self):
        """Initialize settings from the rc file and CLI.

        Default values are loaded from the resource control
        file (processPOPTCrc.py) found in $OLCAO_RC or in
        the current working directory. These are then
        overridden by any explicit command-line arguments.
        """

        # Read default variables from the resource control
        # file.
        rc_dir = os.getenv('OLCAO_RC')
        if not rc_dir:
            sys.exit(
                "Error: $OLCAO_RC is not set. "
                "See instructions."
            )
        sys.path.insert(1, rc_dir)
        from processPOPTCrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the rc
        # file.
        self.reconcile(args)

        # Record the command line parameters that were used.
        self.record_clp()

    def assign_rc_defaults(self, default_rc):
        """Copy rc-file defaults into instance attributes.

        Parameters
        ----------
        default_rc : dict
            Dictionary returned by
            processPOPTCrc.parameters_and_defaults().
        """

        # Spin channel: 1 = up/total, 2 = down.
        self.spin = default_rc["spin"]

    def parse_command_line(self):
        """Build the argument parser and parse sys.argv.

        Returns
        -------
        argparse.Namespace
            Parsed arguments.
        """

        prog_name = "processPOPTC.py"

        description_text = """\
Partial Optical Properties (POPTC) post-processor.

Reads the raw partial Epsilon 2 data file produced by the
OLCAO poptc calculation (fort.250 or fort.251) and processes
each partial pair through OLCAOkkc to produce the
corresponding partial Epsilon 1, ELF, refractive index,
extinction coefficient, Epsilon 1 imaginary, reflectivity,
and absorption coefficient spectra.

Each partial epsilon2 is extracted via makePDOS, run through
OLCAOkkc with the appropriate scaling factor, and the results
are assembled into raw output files for downstream plotting
or analysis.
"""

        epilog_text = """\
Please contact Paul Rulis regarding questions.
Defaults are given in ./processPOPTCrc.py or
$OLCAO_RC/processPOPTCrc.py.
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
        """Register all CLI arguments with the parser.

        Parameters
        ----------
        parser : argparse.ArgumentParser
            The parser to which arguments are added.
        """

        # The -s option selects the spin channel.
        # In the Perl original this was a positional argument
        # ($ARGV[0]). Here we use a named option for clarity.
        # Spin 1 = up or total (non-spin-polarized).
        # Spin 2 = down.
        parser.add_argument(
            '-s', '--spin', dest='spin',
            type=int, default=self.spin,
            choices=[1, 2],
            help=(
                "Spin channel to process. "
                "1 = spin up or total (default). "
                "2 = spin down."
            ),
        )

    def reconcile(self, args):
        """Merge parsed CLI arguments into settings.

        Parameters
        ----------
        args : argparse.Namespace
            Parsed command-line arguments.
        """

        self.spin = args.spin

    def record_clp(self):
        """Append the command line used to the 'command' file.

        This mirrors the Perl original's behaviour of logging
        every invocation for reproducibility.
        """

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


# ============================================================
# POPTCData: data container and processing engine
# ============================================================

class POPTCData():
    """Container for partial optical properties data and the
    processing logic that drives the KKC pipeline.

    This class manages the full processPOPTC workflow:
      1. Read metadata from the raw eps2 file (style, number
         of sequences, energy scale, element pair labels).
      2. Read KKC scaling factors from the control file.
      3. For each element pair (sequence):
         a. Create a makePDOS control file.
         b. Run makePDOS to extract the pair's eps2 data.
         c. Run OLCAOkkc to produce derived spectra.
         d. Append the derived spectra into raw output files.
      4. Clean up intermediate files.

    Attributes
    ----------
    settings : ScriptSettings
        The user-configurable settings object.
    style : int
        The style flag from the raw eps2 file header.
    num_sequences : int
        Number of partial eps2 units (element pairs) in
        the raw file.
    num_energy_points : int
        Number of energy values on the shared energy scale.
    energy_points : list of float
        The energy scale values. Index 0 is unused; the
        physical data starts at index 1 (1-based indexing,
        matching the Perl original).
    element1 : list of str
        Element 1 name for each sequence. Index 0 is unused.
    element2 : list of str
        Element 2 name for each sequence. Index 0 is unused.
    type1 : list of str
        Type number line for element 1 in each sequence.
        Index 0 is unused.
    type2 : list of str
        Type number line for element 2 in each sequence.
        Index 0 is unused.
    pair_kkc_factor : list of float
        KKC scaling factor for each element pair. These
        factors ensure the sum of decomposed spectra equals
        the original total spectrum. Index 0 is unused.

    File name attributes (set during initialization):
        eps2_raw, kkc_control, pdos_control,
        eps2_temp, e1e2elf_temp, eps1_temp, elf_temp,
        refract_temp, extnct_temp, eps1i_temp,
        reflect_temp, absorp_temp,
        e1e2elf_raw, eps1_raw, elf_raw,
        refract_raw, extnct_raw, eps1i_raw,
        reflect_raw, absorp_raw
    """

    def __init__(self, settings):
        """Set up file names and read all metadata.

        Parameters
        ----------
        settings : ScriptSettings
            The user-configurable settings object. The spin
            channel determines which input/output fort files
            are used.
        """

        self.settings = settings

        # Initialize file names based on the spin channel.
        self._init_file_names()

        # Initialize data arrays (1-based indexing with
        # placeholder at index 0, matching the Perl original).
        self.energy_points = [None]  # Index 0 unused.
        self.element1 = [""]         # Index 0 unused.
        self.element2 = [""]         # Index 0 unused.
        self.type1 = [""]            # Index 0 unused.
        self.type2 = [""]            # Index 0 unused.
        self.pair_kkc_factor = [None]  # Index 0 unused.

        # Gather metadata and labels from the raw eps2 file.
        self.get_eps2_raw_metadata()

        # Get the KKC control data (scaling factors for each
        # element pair).
        self.get_kkc_control_data()

    def _init_file_names(self):
        """Initialize all file names based on spin channel.

        The OLCAO code uses Fortran unit numbers as file names
        (fort.NNN). The specific numbers depend on whether we
        are processing spin 1 (up/total) or spin 2 (down).
        """

        spin = self.settings.spin

        # Input raw file from OLCAO poptc calculation.
        if spin == 1:
            self.eps2_raw = "fort.250"
        else:
            self.eps2_raw = "fort.251"

        # Data from OLCAO poptc for running poptc kkc.
        self.kkc_control = "fort.209"

        # Intermediate files (created and deleted per
        # sequence iteration).

        # Control file to use makePDOS on poptc output.
        self.pdos_control = "fort.700"

        # Temp input to kkc (the extracted eps2 for one
        # pair).
        self.eps2_temp = "fort.450"

        # Temp output files from kkc (one pair at a time).
        self.e1e2elf_temp = "fort.500"
        self.eps1_temp = "fort.510"
        self.elf_temp = "fort.520"
        self.refract_temp = "fort.530"
        self.extnct_temp = "fort.540"
        self.eps1i_temp = "fort.550"
        self.reflect_temp = "fort.560"
        self.absorp_temp = "fort.570"

        # Raw output files (assemblages of temp files).
        # These accumulate all partial data across all
        # sequences.
        if spin == 1:
            self.e1e2elf_raw = "fort.300"
            self.eps1_raw = "fort.310"
            self.elf_raw = "fort.320"
            self.refract_raw = "fort.330"
            self.extnct_raw = "fort.340"
            self.eps1i_raw = "fort.350"
            self.reflect_raw = "fort.360"
            self.absorp_raw = "fort.370"
        else:
            self.e1e2elf_raw = "fort.301"
            self.eps1_raw = "fort.311"
            self.elf_raw = "fort.321"
            self.refract_raw = "fort.331"
            self.extnct_raw = "fort.341"
            self.eps1i_raw = "fort.351"
            self.reflect_raw = "fort.361"
            self.absorp_raw = "fort.371"

    def get_eps2_raw_metadata(self):
        """Read metadata from the raw eps2 file.

        The raw eps2 file (fort.250 or fort.251) is produced
        by the OLCAO poptc calculation. Its structure is:

          Line 1: STYLE <style_number>
          Line 2: NUM_UNITS <num_sequences>
          Line 3: NUM_POINTS <num_energy_points>
          Next num_energy_points lines: one energy value each
          Then, for each sequence (element pair):
            ELEMENT_1_NAME <name>
            ELEMENT_2_NAME <name>
            TYPE_1 <number>
            TYPE_2 <number>
            (followed by the eps2 data for that pair)

        This subroutine reads only the metadata (style,
        counts, energy scale, and element/type labels). It
        does not read the actual eps2 data values.
        """

        with open(self.eps2_raw, "r") as f:

            # Get the file style (1, 2, etc.)
            values = prep_line(f)
            self.style = values[-1]

            # Get the number of POPTC units in this data
            # file.
            values = prep_line(f)
            self.num_sequences = int(values[-1])

            # Get the number of energy points for each unit.
            values = prep_line(f)
            self.num_energy_points = int(values[-1])

            # Get the single set of energy values that apply
            # to each unit.
            for i in range(1, self.num_energy_points + 1):
                values = prep_line(f)
                self.energy_points.append(values[0])

            # Get a list of the element name pairs.
            # We scan through the remainder of the file
            # looking for lines that contain
            # "ELEMENT_1_NAME". When found, the next three
            # lines provide the element 2 name, type 1
            # number, and type 2 number.
            j = 1
            for line in f:
                if "ELEMENT_1_NAME" in line:
                    # Get element one name (this line).
                    self.element1.append(line.strip())

                    # Get element two name.
                    line2 = f.readline().strip()
                    self.element2.append(line2)

                    # Get the type numbers.
                    # Get type 1 number.
                    line3 = f.readline().strip()
                    self.type1.append(line3)

                    # Get type 2 number.
                    line4 = f.readline().strip()
                    self.type2.append(line4)

                    j += 1

    def get_kkc_control_data(self):
        """Read KKC scaling factors from the control file.

        The KKC control file (fort.209) is produced by the
        OLCAO poptc calculation. It contains one line per
        element pair, with each line giving the scaling factor
        that ensures the sum of decomposed spectra equals the
        original total spectrum.

        The factor for each pair is the last value on its
        line.
        """

        with open(self.kkc_control, "r") as f:
            # Get the scaling value for each element pair.
            for i in range(1, self.num_sequences + 1):
                values = prep_line(f)
                self.pair_kkc_factor.append(values[-1])

    def write_header(self, file_handle):
        """Write the standard raw file header.

        Each raw output file starts with a header that
        identifies the style, the number of sequences (element
        pairs), the number of energy points, and then lists
        all the energy scale values. This header is identical
        across all eight raw output file types.

        Parameters
        ----------
        file_handle : file object
            An open file to write the header to.
        """

        # Print the real header.
        file_handle.write(
            f"STYLE   {self.style} \n"
        )
        file_handle.write(
            f"NUM_UNITS  {self.num_sequences} \n"
        )
        file_handle.write(
            f"NUM_POINTS   {self.num_energy_points} \n"
        )

        # Print the energy scale.
        for i in range(1, self.num_energy_points + 1):
            file_handle.write(
                f"{float(self.energy_points[i]):10.4f}\n"
            )

    def create_pdos_control(self, seq_num):
        """Create a control file for makePDOS.

        For each sequence (element pair), we need to extract
        the eps2 data for that pair from the raw file. This
        is done by running makePDOS with a control file that
        specifies which sequence number to extract and how to
        label the columns. The control file requests four
        columns: TOTAL, x, y, z.

        Parameters
        ----------
        seq_num : int
            The 1-based sequence number of the current
            element pair being processed.
        """

        with open(self.pdos_control, "w") as f:
            f.write(
                "LABEL : SEQUENCE_NUM : COL_LABELS\n"
            )
            f.write(f"{seq_num} T : {seq_num} : TOTAL\n")
            f.write(f"{seq_num} x : {seq_num} : x\n")
            f.write(f"{seq_num} y : {seq_num} : y\n")
            f.write(f"{seq_num} z : {seq_num} : z\n")

    def call_make_pdos(self):
        """Run makePDOS to extract one pair's eps2 data.

        Calls makePDOS using the control file (fort.700),
        the raw data file (fort.250 or fort.251), and a new
        output file (fort.450) that will contain data that
        can be plotted. The fort.450 file is effectively an
        eps2 file that can be used directly by the OLCAOkkc
        program.

        After makePDOS completes, the control file is removed
        because it is regenerated for each sequence.
        """

        olcao_bin = os.getenv('OLCAO_BIN', '')

        # Construct the makePDOS command. The -i flag
        # provides the control file, -f provides the raw
        # data file, and -o provides the output file.
        make_pdos_cmd = os.path.join(olcao_bin, "makePDOS")
        subprocess.run(
            [
                make_pdos_cmd,
                "-i", self.pdos_control,
                "-f", self.eps2_raw,
                "-o", self.eps2_temp,
            ],
            check=True,
        )

        # Remove the control file (it will be recreated
        # for the next sequence).
        if os.path.exists(self.pdos_control):
            os.remove(self.pdos_control)

    def call_olcao_kkc(self, seq_num):
        """Run OLCAOkkc on the current pair's eps2 data.

        OLCAOkkc performs the Kramers-Kronig conversion to
        produce derived optical spectra (eps1, ELF,
        refractive index, extinction coefficient, eps1i,
        reflectivity, absorption coefficient) from the eps2
        input.

        The program requires the following arguments:
          1. Number of lines in the eps2 temp file
          2. A flag (1)
          3. Another flag (1)
          4. The KKC scaling factor for this element pair

        Parameters
        ----------
        seq_num : int
            The 1-based sequence number of the current
            element pair being processed.

        Note
        ----
        Currently not set up for spin 2 (see OLCAO source
        line 1376 in the original Perl comment). The spin 2
        pathway exists for file naming but the OLCAOkkc call
        uses the same flags regardless of spin.
        """

        # Count the number of lines in the eps2 temp file.
        with open(self.eps2_temp, "r") as f:
            optc_lines = sum(1 for _ in f)

        olcao_bin = os.getenv('OLCAO_BIN', '')

        # Call OLCAOkkc.
        kkc_cmd = os.path.join(olcao_bin, "OLCAOkkc")
        subprocess.run(
            [
                kkc_cmd,
                str(optc_lines),
                "1",
                "1",
                str(self.pair_kkc_factor[seq_num]),
            ],
            check=True,
        )

    def copy_data(self, file_handle, temp_file_name,
                  num_col_labels, label, seq_num):
        """Copy data from a temp file into a raw output file.

        After OLCAOkkc produces its output for one element
        pair, this subroutine appends that pair's data into
        the corresponding accumulated raw output file. It
        writes the sequence header (element names, type
        numbers, column labels) and then copies all data
        lines, stripping the first two columns (which are
        typically an index and the energy value already
        captured in the file header).

        Parameters
        ----------
        file_handle : file object
            An open file for the accumulated raw output.
        temp_file_name : str
            Path to the temp file produced by OLCAOkkc for
            this spectral type.
        num_col_labels : int
            Number of column labels (3 for eps1/eps2/ELF,
            4 for TOTAL/x/y/z).
        label : str
            Space-separated column label string
            (e.g. "Epsilon1 Epsilon2 ELF" or
            "TOTAL x y z").
        seq_num : int
            The 1-based sequence number of the current
            element pair being processed.
        """

        # Write the sequence header to the raw output file.
        file_handle.write(
            f"SEQUENCE_NUM {seq_num}\n"
        )
        file_handle.write(
            f"{self.element1[seq_num]}\n"
        )
        file_handle.write(
            f"{self.element2[seq_num]}\n"
        )
        file_handle.write(
            f"{self.type1[seq_num]}\n"
        )
        file_handle.write(
            f"{self.type2[seq_num]}\n"
        )
        file_handle.write(
            f"COL_LABELS {num_col_labels}\n"
        )
        file_handle.write(f"{label}\n")

        # Open the temp file, skip the header line, and
        # copy all data lines to the raw output file.
        with open(temp_file_name, "r") as temp_in:
            # Read past the header.
            temp_in.readline()

            # Read all actual data lines in and copy to
            # the raw file. For each line, strip the first
            # two whitespace-separated tokens (typically an
            # index and the energy value) and write the
            # remaining values.
            for line in temp_in:
                values = line.split()
                # Remove the first two columns (index and
                # energy), keeping only the spectral data.
                if len(values) > 2:
                    values = values[2:]
                file_handle.write(
                    " ".join(values) + "\n"
                )

        # Remove the temp file after its data has been
        # copied.
        if os.path.exists(temp_file_name):
            os.remove(temp_file_name)

    def process_all_sequences(self):
        """Run the full processing pipeline.

        This is the main driver that loops over all element
        pair sequences and orchestrates the pipeline:
          1. Open all eight raw output files for writing.
          2. Write the standard header to each.
          3. For each sequence:
             a. Create the makePDOS control file.
             b. Run makePDOS to extract the pair's eps2.
             c. Run OLCAOkkc to produce derived spectra.
             d. Copy each derived spectrum into its raw
                output file.
          4. Close all output files.
          5. Clean up any remaining intermediate files.
        """

        # Open all eight raw output files for writing.
        # Each file accumulates one spectral type across
        # all element pair sequences.
        raw_files = {
            "e1e2elf": open(self.e1e2elf_raw, "w"),
            "eps1": open(self.eps1_raw, "w"),
            "elf": open(self.elf_raw, "w"),
            "refract": open(self.refract_raw, "w"),
            "extnct": open(self.extnct_raw, "w"),
            "eps1i": open(self.eps1i_raw, "w"),
            "reflect": open(self.reflect_raw, "w"),
            "absorp": open(self.absorp_raw, "w"),
        }

        try:
            # Write the standard header (style, num_units,
            # num_points, energy scale) to each raw output
            # file.
            for fh in raw_files.values():
                self.write_header(fh)

            # Begin the main processing loop. For each
            # element pair sequence, extract the eps2 data,
            # run the KKC, and accumulate the results.
            for i in range(1, self.num_sequences + 1):
                # Create a control file for makePDOS for
                # the current sequence.
                self.create_pdos_control(i)

                # Call makePDOS on the current segment to
                # extract this pair's eps2 into the temp
                # file (fort.450).
                self.call_make_pdos()

                # Run OLCAOkkc on the current segment to
                # produce all derived spectra from the
                # extracted eps2.
                self.call_olcao_kkc(i)

                # Copy data from each temp output file
                # into the corresponding raw output file
                # and delete the temp files.
                #
                # The e1e2elf file has 3 columns:
                #   Epsilon1, Epsilon2, ELF
                self.copy_data(
                    raw_files["e1e2elf"],
                    self.e1e2elf_temp,
                    3, "Epsilon1 Epsilon2 ELF", i,
                )

                # The remaining files each have 4 columns:
                #   TOTAL, x, y, z
                for key, temp in [
                    ("eps1", self.eps1_temp),
                    ("elf", self.elf_temp),
                    ("refract", self.refract_temp),
                    ("extnct", self.extnct_temp),
                    ("eps1i", self.eps1i_temp),
                    ("reflect", self.reflect_temp),
                    ("absorp", self.absorp_temp),
                ]:
                    self.copy_data(
                        raw_files[key], temp,
                        4, "TOTAL x y z", i,
                    )

        finally:
            # Close all raw output files.
            for fh in raw_files.values():
                fh.close()

        # Clean up any remaining intermediate files.
        # These should already be deleted by copy_data,
        # but this ensures nothing is left behind.
        temp_files = [
            self.eps2_temp,
            self.e1e2elf_temp,
            self.eps1_temp,
            self.elf_temp,
            self.refract_temp,
            self.extnct_temp,
            self.eps1i_temp,
            self.reflect_temp,
            self.absorp_temp,
        ]
        for temp_file in temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)


# ============================================================
# Main entry point
# ============================================================

def main():
    """Main entry point for processPOPTC.

    The workflow is:
      1. Parse settings from the rc file and command line.
      2. Initialize the POPTCData object (reads metadata
         from the raw eps2 file and KKC control data).
      3. Process all element pair sequences through the
         makePDOS -> OLCAOkkc -> accumulate pipeline.
    """

    # Get script settings from a combination of the resource
    # control file and parameters given by the user on the
    # command line.
    settings = ScriptSettings()

    # Initialize the POPTC data object. This reads all
    # metadata from the raw eps2 file and the KKC control
    # data.
    poptc_data = POPTCData(settings)

    # Process all element pair sequences through the
    # pipeline.
    poptc_data.process_all_sequences()


if __name__ == '__main__':
    # Everything before this point was a subroutine
    # definition or a request to import information from
    # external modules. Only now do we actually start
    # running the program. The purpose of this is to allow
    # another python program to import *this* script and
    # call its functions internally.
    main()
