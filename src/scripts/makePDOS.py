#!/usr/bin/env python3

"""makePDOS — Partial Density of States post-processor.

PROGRAM:  makePDOS.py
PURPOSE:  This program will read a raw PDOS data file for the
   purpose of collecting the various data components and
   combining them according to the requests of the command line
   and a user created control file.

AUTHOR:  Paul Rulis (Perl original); Python port by Claude.
LAST MODIFIED:  Mar. 20, 2026

USAGE:
   makePDOS.py [-i CONTROL_FILE] [-o OUTPUT_FILE]
               [-f RAW_DATA_FILE] [-xanes] [-negToZero]
               [-h | --help]

The -i option is used to provide a control file that defines
   how to collect and organize the data from the raw pdos file.
   The format of the control file is given below.

The -o option is used to override the default output file name.
   The default is PDOS.plot.

The -f option is used to override the default raw data file
   name. The default is gs_dos-fb.p.raw.

The -xanes option will alter the default control file (produced
   internally) to have the curves decomposed into s+d and p
   states.

The -negToZero option will make any negative numbers equal to
   zero in the final output file.

To understand the format of the control file it is necessary to
   first understand the nature of the raw input file. The raw
   PDOS input file has a header followed by a series of data
   units.

   1) The header: The first line of the header lists the number
      of data units in the file. The second line of the header
      lists the number of data points in the energy scale for
      each unit. The last part of the header is a list of the
      data points of the energy scale. This is typically a few
      thousand points with each point on a new line.

   2) The data units: Each data unit is divided into two parts.
      The data unit header, and the data unit data.

   2a) The data unit header consists of a series of filter
       name / filter value pairs. The first filter name is
       always "SEQUENCE_NUM" and the last filter name is always
       "COL_LABELS". The meaning of the first filter value
       depends on whether the data was printed according to the
       atom or atomic type. So, it is either the atom number or
       the type number. The last filter value is actually a long
       string of space separated names that correspond to the
       data columns that were printed for this data unit. In the
       case of raw data printed according to atom, there is only
       one data column called "TOTAL". In the case of raw data
       printed according to atomic types there are many data
       columns with names such as "TOTAL 2s 3s 2p 3p 3d".

   2b) The data unit data consists of columns of data
       corresponding to the energy scale from the data file's
       header. The number of columns for each data type depends
       on the type of print out that was requested. In the case
       where the data was printed according to the atoms, there
       is only one column of data per data unit. This column is
       the total PDOS for that atom. In the case where the data
       was printed according to the atomic types, there are many
       columns of data depending on the valence orbitals that
       were included in this type. There is one column for the
       total PDOS of this type, and there is one column for each
       valence orbital of this type.

The format of the control file can be understood as follows:
   1) The first line is used to define the format of the
      remaining lines. It consists of a series of colon
      separated character strings called tags. The first tag is
      always "LABEL" and the last tag is always "COL_LABELS".
      The tags in between can use any name from the lists of
      available filter names present in the raw data file (e.g.
      SEQUENCE_NUM or ELEMENT_NAME). This first line is called
      the "filter definition line".

   2) The remaining lines are also colon separated lists of
      strings, but they will follow the format defined in the
      filter definition line and include the filter values
      instead of the filter names. Each of these lines will
      correspond to one column of output data and thus one
      curve when plotted. An important thing to understand is
      that the filter values can be listed as a set of values
      within the colon separated list. In this way a number of
      data units can potentially match the curve you are
      defining. Some examples are in order.

   2a) Example filter definition line:
       LABEL : ELEMENT_NAME : SPECIES_ID : TYPE_ID : COL_LABELS
       This line indicates that each of the remaining lines
       should contain a set of values defining a label for the
       curve, which elements should contribute to the curve,
       which species of those elements should contribute to the
       curve, and which types of those species should contribute
       to the curve. The last filter name in this line indicates
       that the available data columns should be filtered now.

   2b) Example curve definition matching the above filter
       definition line:
       TOTAL : all : all : all : TOTAL
       This gives the TDOS. All elements, species and types
       will match, and only the total DOS data from the
       available columns for each unit will be included.

   2c) Example curve definition matching the above filter
       definition line:
       Si2_2 : Si : 2 : 2 : TOTAL
       This gives the total DOS for silicon atoms of species 2
       and type 2.

   2d) Example curve definition matching the above filter
       definition line:
       Si_s : Si : all : all : s
       This gives the combined s orbitals for all silicon.

   2e) Example curve definition matching the above filter
       definition line:
       Si_2_3_4 : Si : 2 3 4 : all : TOTAL
       This gives the combined DOS for silicon of species 2, 3,
       and 4.
"""

import argparse as ap
import math
import os
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
        The delimiter string (or regex-style literal) to split
        on. None means split on whitespace.

    Returns
    -------
    list of str
        The tokens produced by the split, with a leading empty
        string removed if present.
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
    """Holds all user-configurable settings for makePDOS.

    The instance variables of this object are the user settings
    that control the program. The variable values are pulled
    from a list that is created within a resource control file
    and that are then reconciled with command line parameters.
    """

    def __init__(self):
        """Initialize settings from the rc file and CLI.

        Default values are loaded from the resource control
        file (makePDOSrc.py) found in $OLCAO_RC or in the
        current working directory. These are then overridden
        by any explicit command-line arguments.
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
        from makePDOSrc import parameters_and_defaults
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
            makePDOSrc.parameters_and_defaults().
        """

        # File name defaults.
        self.control_file = default_rc["control_file"]
        self.raw_data_file = default_rc["raw_data_file"]
        self.output_file = default_rc["output_file"]

        # Behavioural flags.
        self.neg_to_zero = default_rc["neg_to_zero"]
        self.xanes = default_rc["xanes"]

    def parse_command_line(self):
        """Build the argument parser and parse sys.argv.

        Returns
        -------
        argparse.Namespace
            Parsed arguments.
        """

        prog_name = "makePDOS.py"

        description_text = """\
Partial Density of States (PDOS) post-processor.

Reads a raw PDOS data file produced by the OLCAO DOS program
and combines the various data components according to filters
specified in a control file (or using built-in defaults).

The output is a multi-column text file suitable for plotting,
with one column per requested PDOS curve and an energy column.
"""

        epilog_text = """\
Please contact Paul Rulis regarding questions.
Defaults are given in ./makePDOSrc.py or
$OLCAO_RC/makePDOSrc.py.

See the module docstring (makePDOS.py -h) or the source code
for a detailed description of the control file format.
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

        # The -i option is used to provide a control file that
        # defines how to collect and organize the data from the
        # raw PDOS file. The format of the control file is
        # described in the module docstring.
        parser.add_argument(
            '-i', dest='control_file',
            type=str, default=self.control_file,
            help=(
                "Path to the control file that defines how "
                "to collect and organize the data from the "
                "raw PDOS file. Default: use internal "
                "defaults (no control file)."
            ),
        )

        # The -f option is used to override the default raw
        # data file name. The default is gs_dos-fb.p.raw.
        parser.add_argument(
            '-f', dest='raw_data_file',
            type=str, default=self.raw_data_file,
            help=(
                "Path to the raw PDOS data file. "
                f"Default: {self.raw_data_file}"
            ),
        )

        # The -o option is used to override the default output
        # file name. The default is PDOS.plot.
        parser.add_argument(
            '-o', dest='output_file',
            type=str, default=self.output_file,
            help=(
                "Path to the output file. "
                f"Default: {self.output_file}"
            ),
        )

        # The -xanes option will alter the default control file
        # (produced internally) to have the curves decomposed
        # into s+d and p states.
        parser.add_argument(
            '-xanes', dest='xanes',
            action='store_true', default=self.xanes,
            help=(
                "Use XANES-style orbital decomposition "
                "(s+d and p) for the default control. "
                "Only meaningful when no -i control file "
                "is given."
            ),
        )

        # The -negToZero option will make any negative numbers
        # equal to zero in the final output file.
        parser.add_argument(
            '-negToZero', dest='neg_to_zero',
            action='store_true', default=self.neg_to_zero,
            help=(
                "Clamp negative values to zero in the "
                "output file."
            ),
        )

    def reconcile(self, args):
        """Merge parsed CLI arguments into settings.

        Parameters
        ----------
        args : argparse.Namespace
            Parsed command-line arguments.
        """

        self.control_file = args.control_file
        self.raw_data_file = args.raw_data_file
        self.output_file = args.output_file
        self.xanes = args.xanes
        self.neg_to_zero = args.neg_to_zero

    def record_clp(self):
        """Append the command line used to the 'command' file.

        This mirrors the Perl original's behaviour of logging
        every invocation for reproducibility.
        """

        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime("%b. %d, %Y: %H:%M:%S")
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write("Cmnd:")
            for argument in sys.argv:
                cmd.write(f" {argument}")
            cmd.write("\n\n")


# ============================================================
# PDOSData: raw data, control filters, and accumulation
# ============================================================

class PDOSData():
    """Container for raw PDOS data, filter definitions, and
    accumulated output columns.

    This class manages three stages of PDOS processing:
      1. Reading the header (energy scale) from the raw file.
      2. Establishing filter/column definitions (from a
         control file or from internal defaults).
      3. Reading data units, matching them against the
         filters, and accumulating the results.

    Attributes
    ----------
    style : int
        The style flag from the raw PDOS file header (1 =
        per-type output, 2 = pair output, etc.).
    num_units : int
        Number of data units in the raw PDOS file.
    num_energy_points : int
        Number of energy values on the shared energy scale.
    energy_points : list of float
        The energy scale values. Index 0 is unused; the
        values start at index 1 to match the 1-based Perl
        indexing convention.
    format_names : list of str
        The filter format names parsed from the first line
        of the control file (e.g. ["", "LABEL",
        "ELEMENT_NAME", ..., "COL_LABELS"]). Index 0 is
        unused.
    num_col_defs : int
        Number of output column definitions (i.e. the
        number of curves that will appear in the output).
    col_defs : list
        A 3-D structure:
        col_defs[col][fmt_idx] = list of acceptable filter
        term values for output column *col* at format
        position *fmt_idx*. Indices start at 1.
    pdos_col_data : list
        A 2-D structure:
        pdos_col_data[col][energy_idx] = accumulated PDOS
        value for output column *col* at energy index
        *energy_idx*. Indices start at 1.
    using_default_control : bool
        True when the internal default control is used
        (i.e. no explicit control file was provided).
    """

    def __init__(self, settings):
        """Initialize PDOS processing.

        Parameters
        ----------
        settings : ScriptSettings
            The reconciled user settings.
        """

        self.settings = settings

        # These will be populated by init_pdos().
        self.style = 0
        self.num_units = 0
        self.num_energy_points = 0
        self.energy_points = []

        # These will be populated by read_control() or
        # default_control().
        self.format_names = []
        self.num_col_defs = 0
        self.col_defs = []
        self.using_default_control = True

        # This will be populated by read_data().
        self.pdos_col_data = []

    # --------------------------------------------------------
    # init_pdos: read the header of the raw PDOS file
    # --------------------------------------------------------

    def init_pdos(self):
        """Read initialization data from the raw PDOS file.

        Opens the raw PDOS data file and reads the header
        section, which contains:
          - The file style flag (1, 2, etc.)
          - The number of data units in the file
          - The number of energy points in the shared scale
          - The energy scale values themselves

        After this method returns, self.style,
        self.num_units, self.num_energy_points, and
        self.energy_points are populated.
        """

        raw_path = self.settings.raw_data_file

        try:
            self._raw_fh = open(raw_path, 'r')
        except OSError:
            sys.exit(
                f"Cannot open {raw_path} file."
            )

        # Get the file style (1, 2, etc.).
        # The style value is the last token on the line.
        values = prep_line(self._raw_fh)
        self.style = int(values[-1])

        # Get the number of PDOS units in this data file.
        # Again, the count is the last token on the line.
        values = prep_line(self._raw_fh)
        self.num_units = int(values[-1])

        # Get the number of energy points for each unit.
        values = prep_line(self._raw_fh)
        self.num_energy_points = int(values[-1])

        # Get the energy values that apply to each unit.
        # Index 0 is unused (1-based, matching the Perl).
        self.energy_points = [0.0]
        for _ in range(self.num_energy_points):
            values = prep_line(self._raw_fh)
            self.energy_points.append(float(values[0]))

    # --------------------------------------------------------
    # read_control: parse a user-supplied control file
    # --------------------------------------------------------

    def read_control(self):
        """Read and parse a user-supplied control file.

        The control file format is described in the module
        docstring. In brief:
          - The first line is the "filter definition line":
            a colon-separated list of filter category names
            (e.g. "LABEL : ELEMENT_NAME : ... : COL_LABELS").
          - Each subsequent line is a "curve definition":
            a colon-separated list of filter values that
            define one output curve.

        Filter values within a single colon-separated field
        can contain multiple space-separated terms (e.g.
        "2 3 4") meaning *any* of those values will match.
        Numeric ranges (e.g. "2-5") are expanded into the
        individual integers 2, 3, 4, 5.

        After this method returns, self.format_names and
        self.col_defs are populated.
        """

        # Set flag: not using default control.
        self.using_default_control = False

        control_path = self.settings.control_file

        try:
            ctrl_fh = open(control_path, 'r')
        except OSError:
            sys.exit(
                f"Cannot open {control_path} file."
            )

        # -------------------------------------------------
        # Parse the format (filter definition) line.
        # -------------------------------------------------

        # Obtain the format line for the column definitions
        # in the control file. Split on colon.
        values = prep_line(ctrl_fh, splitter=':')

        # Split the format identifiers into separate entries.
        # Index 0 is unused (1-based indexing).
        self.format_names = ['']
        for val in values:
            # Each colon-separated field may contain extra
            # whitespace; take the first non-empty token.
            tokens = val.split()
            if tokens:
                self.format_names.append(tokens[0])

        # -------------------------------------------------
        # Parse each curve (column) definition line.
        # -------------------------------------------------

        self.num_col_defs = 0
        # Pre-fill index 0 so col_defs is 1-based.
        self.col_defs = [None]

        for line in ctrl_fh:
            line = line.strip()
            if not line:
                continue

            # Increment the number of column definitions
            # read in so far.
            self.num_col_defs += 1

            # Read the column definition on this line and
            # store its filter parts in the definition list.
            # Split the line on colon, then strip each part.
            parts = line.split(':')
            parts = [p.strip() for p in parts]

            # Build the filter terms for this column def.
            # col_defs[col_num] is a list (1-based) of lists.
            # col_defs[col_num][fmt_idx] is a list of
            # acceptable values for that format position.
            col_def_entry = [None]  # index 0 unused

            for part in parts:
                # Initialize a list of filter terms for
                # this format position.
                filter_terms = [None]

                # Separate the individual terms.
                term_tokens = part.split()

                for token in term_tokens:
                    if not token:
                        continue

                    # Determine if the filter term is a
                    # numeric range (e.g. "2-5") and expand
                    # it accordingly.
                    if '-' in token:
                        # Try to interpret as a range.
                        range_parts = token.split('-')
                        try:
                            range_start = int(
                                range_parts[0]
                            )
                            range_end = int(
                                range_parts[1]
                            )
                            for k in range(
                                range_start,
                                range_end + 1,
                            ):
                                filter_terms.append(str(k))
                            continue
                        except (ValueError, IndexError):
                            # Not a numeric range; treat as
                            # a literal string (e.g. an
                            # element name or orbital label
                            # that happens to contain '-').
                            pass

                    # Store the individual filter term.
                    filter_terms.append(token)

                col_def_entry.append(filter_terms)

            self.col_defs.append(col_def_entry)

        ctrl_fh.close()

    # --------------------------------------------------------
    # default_control: set up internal default filters
    # --------------------------------------------------------

    def default_control(self):
        """Define default filter and column definitions.

        When no control file is given, this method creates
        a default set of column definitions. The behaviour
        depends on whether XANES decomposition was
        requested:

        Without -xanes:
          One output column per data unit, each showing the
          TOTAL PDOS for that unit.

        With -xanes:
          Two output columns per data unit: one for s+d
          orbitals and one for p orbitals. The first
          num_units columns are s+d; the next num_units
          columns are p.

        The column names (labels) are left blank here and
        are filled in later when the data is read and the
        element/species/type information becomes available.

        After this method returns, self.format_names,
        self.num_col_defs, and self.col_defs are populated.
        """

        # Set flag: using default control.
        self.using_default_control = True

        # Define the format for default columns:
        # format_names[1] = "NAME" (the label, to be filled
        #   in later)
        # format_names[2] = "SEQUENCE_NUM" (which data unit)
        # format_names[3] = "COL_LABELS" (which data columns)
        self.format_names = ['', 'NAME', 'SEQUENCE_NUM',
                             'COL_LABELS']

        # Use the number of units as the default number of
        # columns except for the case of XANES decomposition
        # where we double it. 1 curve for s+d and 1 curve for
        # p states.
        if not self.settings.xanes:
            self.num_col_defs = self.num_units
        else:
            self.num_col_defs = self.num_units * 2

        # Pre-fill index 0 so col_defs is 1-based.
        self.col_defs = [None]

        for i in range(1, self.num_col_defs + 1):
            # Each column definition has 3 format positions
            # (matching format_names above):
            #   [1] = label (to be filled in later)
            #   [2] = sequence number filter terms
            #   [3] = column label filter terms
            col_def_entry = [
                None,  # index 0 unused
                [],    # [1] NAME — filled in during read
                [],    # [2] SEQUENCE_NUM
                [],    # [3] COL_LABELS
            ]

            if not self.settings.xanes:
                # Use the sequence number (i) and request
                # the TOTAL PDOS for this unit.
                col_def_entry[2] = [None, str(i)]
                col_def_entry[3] = [None, 'TOTAL']
            else:
                # The first num_units columns are for the
                # s+d orbitals and the next num_units columns
                # are for the p orbitals.
                if i <= self.num_units:
                    col_def_entry[2] = [None, str(i)]
                    col_def_entry[3] = [
                        None, 's', 'd'
                    ]
                else:
                    col_def_entry[2] = [
                        None,
                        str(i - self.num_units)
                    ]
                    col_def_entry[3] = [None, 'p']

            self.col_defs.append(col_def_entry)

    # --------------------------------------------------------
    # read_data: read units and accumulate matching results
    # --------------------------------------------------------

    def read_data(self):
        """Read PDOS data units and accumulate results.

        This is the core processing routine. For each data
        unit in the raw PDOS file, it:

        1. Reads the unit's filter-name/filter-value header
           pairs (e.g. SEQUENCE_NUM, ELEMENT_NAME, etc.)
           into a dictionary (unit_hash).

        2. Reads the COL_LABELS line to learn which data
           columns are available in this unit (e.g.
           "TOTAL 2s 3s 2p 3p 3d").

        3. Reads the numerical data columns for this unit
           (one row per energy point, one column per label).

        4. Compares each output column definition against
           the unit's filter values. If the unit matches
           *all* filter terms for a given column definition,
           then each requested data column label is compared
           to the available column labels. Whenever a match
           is found, that column of unit data is accumulated
           into the output column.

        5. When using default control settings, assigns the
           column name (label) at the time of the first
           match, using information from the unit hash (e.g.
           element name, species ID, type ID).

        After this method returns, self.pdos_col_data
        contains the accumulated PDOS values for all output
        columns.
        """

        # Initialize the PDOS data accumulators.
        # pdos_col_data[col][energy_idx] = 0.0
        # Both indices are 1-based.
        self.pdos_col_data = [None]  # index 0 unused
        for _ in range(1, self.num_col_defs + 1):
            col_data = [0.0]  # index 0 unused
            for _ in range(self.num_energy_points):
                col_data.append(0.0)
            self.pdos_col_data.append(col_data)

        # Read each PDOS data unit.
        for unit_idx in range(1, self.num_units + 1):

            # -------------------------------------------------
            # Read the unit header (filter name/value pairs).
            # -------------------------------------------------

            # Repeatedly read a filter name from the PDOS
            # file and store it in the hash for this unit.
            unit_hash = {}
            input_col_labels = ['']  # 1-based

            while True:
                values = prep_line(self._raw_fh)

                # Get the key name.
                current_key = values[0]

                # Get the value. In the case of COL_LABELS
                # we save the values in a list; all other
                # cases are stored in a dictionary.
                if current_key == 'COL_LABELS':
                    num_col_labels = int(values[1])

                    # Column labels are printed 6 per line.
                    if num_col_labels % 6 == 0:
                        num_lines = num_col_labels // 6
                    else:
                        num_lines = (
                            num_col_labels // 6 + 1
                        )

                    # Obtain all the column label values
                    # and then leave the while loop.
                    for _ in range(num_lines):
                        vals = prep_line(self._raw_fh)
                        for v in vals:
                            input_col_labels.append(v)
                    break
                else:
                    unit_hash[current_key] = values[1]

            # -------------------------------------------------
            # Read the data columns for this unit.
            # -------------------------------------------------

            # unit_data[energy_idx] = list of column values
            # (1-based for both indices).
            unit_data = [None]  # index 0 unused
            for energy_idx in range(
                1, self.num_energy_points + 1
            ):
                row = ['']  # 1-based
                for _ in range(num_lines):
                    vals = prep_line(self._raw_fh)
                    for v in vals:
                        row.append(float(v))
                unit_data.append(row)

            # -------------------------------------------------
            # Match this unit against each column definition.
            # -------------------------------------------------

            # Consider each column definition in turn and
            # determine if this data unit matches the filter
            # terms appropriately. When it matches all the
            # filter terms, then each requested data column
            # from the column definition will be processed.
            # This means that each requested column will be
            # compared to the available columns. When an
            # available column in the provided data matches a
            # requested column, then that column of unit data
            # will be accumulated into the data for the output
            # column of that column definition.

            for col_idx in range(
                1, self.num_col_defs + 1
            ):
                # Assume that the data unit matches this
                # definition.
                matches_definition = True

                # Consider all column definition format names
                # except the output label (index 1) and the
                # input data column label (last index). The
                # current unit must match every definition
                # format name to have the chance of being
                # included in the current column's output.
                num_formats = len(self.format_names) - 1

                for fmt_idx in range(2, num_formats):
                    # If the data unit does not match the
                    # definition, abort and move on to the
                    # next column definition.
                    if not matches_definition:
                        break

                    # Take the current format name and use
                    # it as a key to get the value of that
                    # term for this input data unit.
                    fmt_name = self.format_names[fmt_idx]
                    current_term_value = unit_hash.get(
                        fmt_name, ''
                    )

                    # Compare this data unit's term value to
                    # the acceptable values for this column
                    # definition's term.
                    acceptable = (
                        self.col_defs[col_idx][fmt_idx]
                    )

                    # Check if it matches any acceptable
                    # value (or the wildcard "all").
                    matched_term = False
                    for acc_val in acceptable[1:]:
                        if (current_term_value == acc_val
                                or acc_val == 'all'):
                            matched_term = True
                            break

                    if not matched_term:
                        # The term value for this data unit
                        # does not match any acceptable terms
                        # for this column.
                        matches_definition = False

                # If there was no match, move on to the next
                # column definition.
                if not matches_definition:
                    continue

                # If there was a match so far, now we compare
                # the available data columns for this unit
                # with the requested column labels. Any time
                # a column matches a request it will be added
                # to the PDOS for this output column.
                last_fmt = len(self.format_names) - 1
                requested = self.col_defs[col_idx][last_fmt]

                for req_label in requested[1:]:
                    # Compare the requested label to each of
                    # the available columns.
                    for label_idx in range(
                        1, len(input_col_labels)
                    ):
                        # Check for a match. The Perl code
                        # uses index() which checks for a
                        # substring match — we replicate
                        # that with Python's "in" operator.
                        if (req_label in
                                input_col_labels[label_idx]
                                or req_label == 'all'):
                            # Accumulate the data.
                            for m in range(
                                1,
                                self.num_energy_points + 1,
                            ):
                                self.pdos_col_data[
                                    col_idx
                                ][m] += unit_data[m][
                                    label_idx
                                ]

                            # If using the default column
                            # definitions then we need to
                            # give the column a name. Since
                            # the default uses only one data
                            # unit per column we can give
                            # the name now that the match
                            # has been found.
                            if self.using_default_control:
                                self._assign_default_name(
                                    col_idx,
                                    unit_hash,
                                )

    def _assign_default_name(self, col_idx, unit_hash):
        """Assign a column name when using default control.

        When no explicit control file is given, the column
        name is derived from the unit's metadata. The exact
        format depends on the file style and whether XANES
        decomposition is active.

        Parameters
        ----------
        col_idx : int
            The 1-based index of the output column.
        unit_hash : dict
            The filter name/value pairs for the current
            data unit.
        """

        if not self.settings.xanes:
            tail_id = ''
        else:
            if col_idx <= self.num_units:
                tail_id = '_sd'
            else:
                tail_id = '_p'

        if self.style == 1:
            name = (
                unit_hash['ELEMENT_NAME']
                + '_'
                + unit_hash['SPECIES_ID']
                + '_'
                + unit_hash['TYPE_ID']
                + tail_id
            )
        elif self.style == 2:
            name = (
                unit_hash['ELEMENT_1_NAME']
                + '_'
                + unit_hash['ELEMENT_2_NAME']
                + tail_id
            )
        else:
            sys.exit(
                f"Unknown style {self.style}. "
                "Stopping."
            )

        self.col_defs[col_idx][1] = [None, name]

    # --------------------------------------------------------
    # print_results: write the output file
    # --------------------------------------------------------

    def print_results(self):
        """Write the accumulated PDOS data to the output file.

        Produces a whitespace-delimited text file with:
          - A header line of column labels (ENERGY followed
            by one label per output column).
          - One row per energy point, with the energy value
            followed by the accumulated PDOS value for each
            output column.

        If the neg_to_zero setting is active, any negative
        PDOS values are replaced with 0.0 in the output.
        """

        output_path = self.settings.output_file

        try:
            out_fh = open(output_path, 'w')
        except OSError:
            sys.exit(
                f"Cannot open {output_path} file."
            )

        # Put on the header of column labels.
        out_fh.write(f"{'ENERGY':>12s}")
        for col_idx in range(1, self.num_col_defs + 1):
            label = self.col_defs[col_idx][1][1]
            out_fh.write(f"{label:>12s}")
        out_fh.write('\n')

        # Write one row per energy point.
        for energy_idx in range(
            1, self.num_energy_points + 1
        ):
            # Print the energy point for this iteration.
            out_fh.write(
                f"{self.energy_points[energy_idx]:12.4f}"
            )

            for col_idx in range(
                1, self.num_col_defs + 1
            ):
                val = self.pdos_col_data[col_idx][
                    energy_idx
                ]
                if self.settings.neg_to_zero and val < 0:
                    val = 0.0
                out_fh.write(f"{val:12.4f}")

            out_fh.write('\n')

        out_fh.close()


# ============================================================
# Main
# ============================================================

def main():
    """Main entry point for makePDOS processing.

    Execution flow:
      1. Initialize program environment from command line
         parameters and the resource control file.
      2. Read the initialization data (energy scale) from
         the raw PDOS file.
      3. Read the control file if one exists to define the
         filters. If no control file was given, then define
         default filters.
      4. Read the PDOS data file and accumulate the column
         results.
      5. Print the results.
    """

    # Get script settings from a combination of the resource
    # control file and parameters given by the user on the
    # command line.
    settings = ScriptSettings()

    # Create the data container and begin processing.
    pdos = PDOSData(settings)

    # Read the initialization data from the raw PDOS file.
    pdos.init_pdos()

    # Read the control file if one exists to define the
    # filters. If no control file was given, then define
    # default filters.
    if settings.control_file:
        pdos.read_control()
    else:
        pdos.default_control()

    # Read the PDOS data file and accumulate the column
    # results.
    pdos.read_data()

    # Print the results.
    pdos.print_results()


if __name__ == '__main__':
    # Everything before this point was a subroutine
    # definition or a request to import information from
    # external modules. Only now do we actually start
    # running the program. The purpose of this is to allow
    # another python program to import *this* script and
    # call its functions internally.
    main()
