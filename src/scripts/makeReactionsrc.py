#!/usr/bin/env python3

"""makeReactionsrc.py -- Resource control file for makeReactions.py.

This file defines the default parameter values for the makeReactions.py
script.  Users can override these defaults via command line arguments.
A local copy of this file (in the current working directory) takes
precedence over the copy in $OLCAO_RC.

Each parameter is documented with its purpose, valid values, and how
it relates to the physical / algorithmic behaviour of the script.
"""


def parameters_and_defaults():
    """Return a dictionary of default parameters for makeReactions.py.

    Returns
    -------
    dict
        Keys are parameter names; values are the defaults.
    """
    param_dict = {

        # ----- Chain length parameters -----

        # chain_len1 : int
        #   Number of bond-hops from the surface (S) atom of
        #   molecule 1 to include in the reaction template.
        #   A value of 4 captures the S atom, the atom bonded
        #   to S, its nearest neighbours, and the next-nearest
        #   neighbours.  This depth is needed so that bond
        #   angles involving the bonded atom are fully
        #   represented.
        "chain_len1": 4,

        # chain_len2 : int
        #   Same as chain_len1 but for molecule 2.
        "chain_len2": 4,

        # ----- Delete-tree parameters -----

        # del_tree_depth1 : int
        #   Depth of the delete tree for molecule 1.  When 0,
        #   only the S atom itself is deleted upon binding.
        #   When > 0, atoms between the S atom and the target
        #   element/species are also deleted.
        "del_tree_depth1": 0,

        # del_tree_depth2 : int
        #   Same as del_tree_depth1 but for molecule 2.
        "del_tree_depth2": 0,

        # del_tree_target1 : str
        #   Element_species string (e.g. "si_1") that marks
        #   the end of the delete tree for molecule 1.  Only
        #   meaningful when del_tree_depth1 > 0.
        "del_tree_target1": "",

        # del_tree_target2 : str
        #   Same as del_tree_target1 but for molecule 2.
        "del_tree_target2": "",

        # ----- Surface-atom parameters -----

        # s_element : str
        #   Element name (lower-case) that serves as the
        #   "surface" (S) atom.  These are the atoms that
        #   will be removed (or bridged) when two molecules
        #   react.
        "s_element": "h",

        # ss_dist : float
        #   The distance (Angstroms) separating the two S
        #   atoms in the reaction template.
        "ss_dist": 0.5,

        # del_both_s : int  (0 or 1)
        #   1 = delete both S atoms upon reaction.
        #   0 = delete only the S atom from molecule 2;
        #       the S atom from molecule 1 remains and serves
        #       as the bridge between the two trigger atoms.
        "del_both_s": 1,

        # ----- Type-numbering parameters -----

        # type_step : int
        #   Increment applied to atom type numbers upon
        #   bonding.  0 means types remain unchanged
        #   throughout the simulation.
        "type_step": 0,

        # const_types : bool
        #   If True, atom type numbers from precursor
        #   skeleton files are kept constant when molecules
        #   are merged (e.g. C type 1 in mol1 and C type 1
        #   in mol2 remain type 1 in the merged system).
        #   If False, types from molecule 2 are shifted
        #   upward so they are distinct from molecule 1.
        "const_types": False,
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
