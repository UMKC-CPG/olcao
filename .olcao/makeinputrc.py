#!/usr/bin/env python3

"""makeinputrc.py -- Default parameters for makeinput.py.

This resource control file provides a single function that returns a
dictionary of all default parameters used by makeinput.  Users may place
a customized copy of this file in the current working directory or in
the directory pointed to by $OLCAO_RC to override any defaults.
"""

import os


def parameters_and_defaults():
    # Determine queue type from environment (0=bash, 1=pbs, 2=lsf, 3=slurm).
    queue_env = os.getenv("OLCAO_QUEUE", "0")
    try:
        queue_type = int(queue_env)
    except ValueError:
        queue_type = 0

    param_dict = {

        # ---- Executable names ----
        "makekpoints_exec": "makekpoints",
        "contract_exec":    "contract",

        # ---- Potential modification (off by default) ----
        "mod_pot":          0,       # 0=no modification
        "mod_element_name": "",
        "min_mod_term":     0,
        "max_mod_term":     0,
        "num_mod_terms":    0,

        # ---- Output file flags ----
        "pdb":              0,       # 1=produce PDB file
        "cif":              0,       # 1=produce CIF file
        "queue_type":       queue_type,

        # ---- Cutoff criteria ----
        "bf_cutoff":        16,      # Basis function interaction exponent
        "es_cutoff":        16,      # Electrostatic interaction exponent

        # ---- Exchange correlation mesh ----
        "xc_code":          100,     # Wigner interpolation
        "xc_in_weight":     0.5,
        "xc_out_weight":    0.5,
        "xc_in_samp":       0.1,
        "xc_out_samp":      3.5,
        "xc_spacing_samp":  0.8,
        "num_samp_vectors": 100,

        # ---- DOS defaults ----
        "e_delta_dos":      0.01,
        "sigma_dos":        0.5,
        "emin_dos":         -30,
        "emax_dos":         30,
        "detail_code_pdos": 0,       # 0=type, 1=a total, 2=a nl, 3=a nlm
        "iter_tdos":        0,

        # ---- BOND defaults ----
        "max_len_bond":     3.5,
        "emin_bond":        -30,
        "emax_bond":        30,
        "sigma_bond":       0.1,
        "e_delta_bond":     0.01,
        "output_bondq":     0,       # 0=old(type), 1=new(atom)
        "bond_3c":          0,       # 0=no 3-center bond order
        "max_neighbor_bond": 20,

        # ---- PACS (ELNES/XANES) defaults ----
        "e_delta_pacs":     0.01,
        "sigma_pacs":       0.5,
        "onset_slack_pacs": 5.0,
        "energy_window_pacs": 50.0,

        # ---- OPTC defaults ----
        "e_delta_optc":     0.01,
        "e_delta_nlop":     0.01,
        "e_delta_sige":     0.001,
        "sigma_optc":       0.1,
        "sigma_nlop":       0.1,
        "sigma_sige":       0.1,
        "e_trans_optc":     50,
        "e_trans_nlop":     50,
        "e_trans_sige":     0.3,
        "e_cutoff_optc":    45,
        "e_cutoff_nlop":    45,
        "e_cutoff_sige":    5,
        "detail_code_poptc": 0,      # 0=NONE,1=elem,2=a tot,3=e nl,4=e nlm

        # ---- MAIN (SCF) defaults ----
        "num_iter_main":    50,
        "converg_main":     0.0001,
        "spin_main":        0.01,
        "therm_smear_main": 0.0,

        # ---- SYBD / cell name ----
        "sybd_path":        "",
        "cell_name":        "",

        # ---- K-point defaults ----
        "kp_mesh_scf":      [1, 1, 1],
        "kp_mesh_pscf":     [1, 1, 1],
        "kp_weight_sum":    2,
        "kp_shift":         "-1 -1 -1",
        "print_bz":         0,
        "scale_factor":     0.0,

        # ---- Basis/potential substitution (none by default) ----
        "num_basis_subs":   0,
        "num_pot_subs":     0,

        # ---- Grouping methods (all off by default) ----
        # reduce defaults (per-use)
        "reduce_level":     2,
        "reduce_thick":     0.10,
        "reduce_cutoff":    4.0,
        "reduce_op":        "species",
        "reduce_tolerance": 0.05,
        "reduce_selection": 0,
        # target defaults (per-use)
        "target_radius":    3.50,
        "target_zone":      "in",
        "target_op":        "species",
        "target_relation":  "diff",
        # block defaults (per-use)
        "block_zone":       "in",
        "block_op":         "species",
        "block_relation":   "diff",

        # ---- XANES defaults ----
        "xanes":            0,       # 0=off
        "xanes_radius":     3.50,

        # ---- Misc defaults ----
        "state_factor":     2.5,
        "rel":              0,       # 0=non-relativistic
        "do_basis_vis":     0,
        "emu":              0,
        "no_core":          0,       # 0=core orbitals orthogonalized out

        # ---- Slurm defaults ----
        "partition":        "rulisp-lab",
        "account":          "rulisp-lab",
        "time":             "00:60:00",
        "memory":           "10G",
        "cpus":             1,
        "nodes":            1,
    }
    return param_dict


if __name__ == "__main__":
    for k, v in parameters_and_defaults().items():
        print(f"  {k:24s} = {v!r}")
