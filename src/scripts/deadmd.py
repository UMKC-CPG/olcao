#!/usr/bin/env python3
"""
=============================================================================
DEAD_MD: Dynamic Evolutionary Algorithm-Driven Molecular Dynamics
=============================================================================

Version: 1.0
Last Updated: March 2026

Copyright (c) 2026
Mohammed Belhadj Larbi

Author:      Mohammed Belhadj Larbi
Developer:   Mohammed Belhadj Larbi
Institution: University of Missouri-Kansas City (UMKC)

This code was fully designed and implemented by Mohammed Belhadj Larbi.

-----------------------------------------------------------------------------
Description
-----------------------------------------------------------------------------
DEAD-MD is a genetic algorithm framework designed to optimize LAMMPS 
bond/react condensation parameters in order to generate improved amorphous
structural models.

The method couples molecular dynamics simulations with evolutinary search,
allowing model parameters to evolve across generations according to physically
meaningful fitness criteria derived from simulation outputs.

-----------------------------------------------------------------------------
Algorithm Overview
-----------------------------------------------------------------------------
The algorithm operates by:

1. Initialize a population of candidates genomes, where each genome encodes
a set of LAMMPS input parameters.

2. Setting and running LAMMPS condensation simulations for each population member.

3. Parsing simulation outputs to evaluate fitness based on quantities such 
as total energy, density, pair distribution function (PDF), and target
element content percentage.

4. Selecting parents according to the chosen selection strategy.

5. Generating offspring through crossover operations.

6. applying mutation to introduce diversity in the population.

7. Constructing a new generation from surviving individuals and offspring.

8. Iterating the process until the maximum generation limit
is reached or convergence criteria are satisfied.

----------------------------------------------------------------------------
Notes
----------------------------------------------------------------------------
This script contains the complete autonomous implementation of the DEAD-MD
workflow, automatically performing: population initialization, simulation
orchestration, fitness evaluation, parents selection, crossover, mutation,
and generation management.

----------------------------------------------------------------------------
Example deadmd.in
----------------------------------------------------------------------------
The following shows a complete annotated input file.  Every keyword is
case-insensitive.  Lines are parsed by splitting on whitespace; order
within the file does not matter.

    num_cores          4
    gen_max            20
    population_size    10
    elitism_rate       2
    mutation_rate      5
    kill_rate          4
    crossover_type     1
    selection          tournament 3

    target_element     h
    target_element_lower  30.0
    target_element_upper  40.0

    target_density_lower  1.8
    target_density_upper  2.4

    composition_num    2
    c6h19nsi2_1  c-1  50
    c2h6_1       c-1  10

    cell_size_lower    20.0
    cell_size_upper    40.0

    max_speed          5.0

    # ref_energy is the total energy (in whatever units LAMMPS reports)
    # of the best-known or target structure.  It anchors the energy term
    # of the fitness function so that scores are comparable across all
    # generations.  A member whose energy equals ref_energy scores 0 on
    # this term; worse structures score positive; better structures score
    # negative (rewarded).  Use a DFT or force-field minimization result
    # for the composition of interest as a sensible starting value.
    ref_energy         -85432.7

    # weight_energy, weight_density, and weight_element scale the three
    # fitness terms relative to each other.  All three default to 1.0
    # when omitted, which reproduces the original equal-weight behaviour.
    # Raise a weight to make DEAD-MD care more about that quantity; set
    # it to 0.0 to ignore the term entirely.
    weight_energy      1.0
    weight_density     2.0
    weight_element     1.5

    # exp_pdf_file provides the path to an experimental G(r) measured by
    # neutron (or X-ray) diffraction.  The file must contain exactly two
    # rows: the first row is the distance grid in Angstrom (0.01 spacing),
    # the second row is G(r).  When this keyword is present and weight_pdf
    # is greater than zero, a scale-optimised R-factor between the simulated
    # and experimental G(r) is added as a fourth fitness term.  Omitting
    # either keyword disables the PDF term entirely and the script behaves
    # exactly as it did before (fully backward-compatible).
    exp_pdf_file       gr_neutron.dat
    weight_pdf         1.0

    stage              100 nvt 300.0 300.0 100.0 50000


    reactions_num 2  1
    c6h19nsi2_1 c-1 c2h6_1 c-1 0.2 0.8

"""

import os
import sys
import numpy as np
import random
import subprocess


def print_banner():
        print("\n")
        print("==========================================================")
        print(" DEAD-MD: Dynamic Evolutionary Algorithm-Driven Molecular Dynamics")
        print(" Version: 1.0")
        print(" Author:  Mohammed Belhadj Larbi")
        print(" Institution: University of Missouri-Kansas City (UMKC)")
        print("==========================================================")
        print("\n")


def initialize_population():
    """ 
    Reads the input file deadmd.in, and then based on that it will create 
    an initial population of N LAMMPS simulations with their necessary 
    files (data files, input files, reaction templates).

    Args:
        None

    Returns:
        Population, population size, mutation rate, kill rate, elitism rate,
        target element, number of cores, selection method, maximum generation
        number, target_element_lower, target element lower and upper bounds,
        target density lower and upper bounds, cell size lower and upper bounds,
        bonding probabilities lower and upper bounds, crossover type used, max
        molecular speed, simulation box squishing factor, box squishing steps,
        type of ensemble, ensemble initial temperature, ensemble final
        temperature, ensemble temperature rate, simulation number of steps,
        tournament size k, the user-supplied reference total energy (ref_energy)
        used by the fitness function to anchor energy scores across generations,
        and the three fitness term weights (weight_energy, weight_density,
        weight_element) that scale the energy, density, and element-percentage
        contributions to the total fitness score. All three weights default to
        1.0 when the corresponding keyword is omitted from deadmd.in.
        exp_pdf is a numpy array of the experimental G(r) values truncated
        to 0.01-10.00 Angstrom (1000 points matching rpdf's output grid),
        or None when exp_pdf_file is absent or weight_pdf is 0. weight_pdf
        scales the PDF R-factor term and defaults to 0.0 so the term is
        inactive unless the user explicitly enables it in deadmd.in.
    """
    population = []
    molecules = []
    molecule_family = []
    num_molecule = []
    rxns = []
    k = 2 # default value if user does not provide it or choses a selection
          # method that does not require a k value

    # Fitness term weights default to 1.0 (equal weighting). The user may
    # override any or all of them with weight_energy / weight_density /
    # weight_element in deadmd.in without changing the other terms.
    # weight_pdf defaults to 0.0: the PDF term is inactive unless the
    # user supplies both exp_pdf_file and a positive weight_pdf.
    weight_energy  = 1.0
    weight_density = 1.0
    weight_element = 1.0
    weight_pdf     = 0.0
    exp_pdf_file   = None

# Read the deadmd.in input file
    with open("deadmd.in", "r") as input_file:
        input_file_lines = [line.strip() for line in input_file if line.strip()]

        for line_num, line in enumerate(input_file_lines):
            words = line.split()
            if (words[0].lower() == "num_cores"):
                num_cores = int(words[1])
            elif (words[0].lower() == "gen_max"):
                gen_max = int(words[1])
            elif (words[0].lower() == "population_size"):
                population_size = int(words[1])
            elif (words[0].lower() == "elitism_rate"):
                elitism_rate = int(words[1])
            elif (words[0].lower()) == "mutation_rate":
                mutation_rate = float(words[1])/100
            elif (words[0].lower() == "kill_rate"):
                kill_rate = int(words[1])
            elif (words[0].lower() == "crossover_type"):
                crossover_type = int(words[1])
            elif (words[0].lower() == "selection"):
                selection_method = words[1]
                if selection_method == "truncation":
                    selection = 1
                elif selection_method == "rank":
                    selection = 2
                elif selection_method == "tournament":
                    selection = 3
                    k = int(words[2])
            elif (words[0].lower() == "target_element"):
                target_element = words[1].lower()
            elif (words[0].lower() == "target_element_lower"):
                target_element_lower = float(words[1])
            elif (words[0].lower() == "target_element_upper"):
                target_element_upper = float(words[1])
            elif (words[0].lower() == "target_density_lower"):
                target_density_lower = float(words[1])
            elif (words[0].lower()) == "target_density_upper":
                target_density_upper = float(words[1])
            elif (words[0].lower() == "composition_num"):
                composition_num = int(words[1])

                for molecule_line_index in range(line_num + 1, line_num + 1 + 
                        composition_num):
                    parts = input_file_lines[molecule_line_index].split()
                    molecules.append(parts[0])
                    molecule_family.append(parts[1])
                    num_molecule.append(parts[2])
                print(molecules)
                print(molecule_family)
                print(num_molecule)
            elif (words[0].lower() == "cell_size_lower"):
                cell_size_lower = float(words[1])
            elif (words[0].lower() == "cell_size_upper"):
                cell_size_upper = float(words[1])
            elif (words[0].lower() == "max_speed"):
                max_speed = float(words[1])
            elif (words[0].lower() == "stage"):
                squish_step_size = int(words[1])
                ensemble_type = words[2]
                ensemble_t_start = float(words[3])
                ensemble_t_end = float(words[4])
                ensemble_t_rate = float(words[5])
                run_steps = int(words[6])
            elif (words[0].lower() == "ref_energy"):
                ref_energy = float(words[1])
            elif (words[0].lower() == "weight_energy"):
                weight_energy = float(words[1])
            elif (words[0].lower() == "weight_density"):
                weight_density = float(words[1])
            elif (words[0].lower() == "weight_element"):
                weight_element = float(words[1])
            elif (words[0].lower() == "exp_pdf_file"):
                exp_pdf_file = words[1]
            elif (words[0].lower() == "weight_pdf"):
                weight_pdf = float(words[1])
            elif (words[0].lower() == "reactions_num"):
                reactions_num = int(words[1])

                for rxn_line_index in range(line_num + 1, line_num +
                        1 + reactions_num):
                    rxns.append(input_file_lines[rxn_line_index].split())
                bonding_probability_lower = float(rxns[0][4])
                bonding_probability_upper = float(rxns[0][5])


        # Filling the population array. The structure is a 2d array. 
        # Each row is a member (individual or genome).
        # Each column = [target_element_percentage, target_density, 
        #  cell_size, condensation_rate, bonding_probability]
        # The bonding probabilities depend on how many reactions you have, 
        #  if you have 1 reaction than the size of the "population" array is 
        #  (poulation_size, 5), if you have two reactions, then the column size
        #  increases to 6 instead of 5 (population_size, 6) and so on.
        # The sequence of probabilities in the array is in respect to the sequence
        #  of the reaction you list in the deadmd.in or lamps.in
        population = np.empty((population_size, 5))
        for i in range(population_size):
            population[i][0] = random.uniform(target_element_lower, target_element_upper)
            population[i][1] = random.uniform(target_density_lower, target_density_upper)
            population[i][2] = round(random.uniform(cell_size_lower, cell_size_upper), 2)
            population[i][3] = round(random.uniform(0.25, 0.99), 2)
            population[i][4] = random.uniform(bonding_probability_lower, 
                    bonding_probability_upper)
        print(population)

        # Load the experimental PDF if the user supplied a file and a
        # positive weight. The file is two-column whitespace-delimited:
        # column 0 is r (Angstrom), column 1 is G(r). Points are kept
        # up to 10.00 Angstrom to match rpdf's output grid (1000 pts).
        exp_pdf = None
        if exp_pdf_file is not None and weight_pdf > 0.0:
            data    = np.loadtxt(exp_pdf_file)
            r_exp   = data[:, 0]
            g_exp   = data[:, 1]
            mask    = r_exp <= 10.0 + 1e-9
            exp_pdf = g_exp[mask]
            print(f"[DEAD-MD] Experimental PDF loaded: "
                  f"{exp_pdf.size} points (0.01 to 10.00 Angstrom)")

    return (
                population,
                population_size,
                mutation_rate,
                kill_rate,
                elitism_rate, 
                target_element,
                num_cores, 
                selection, 
                gen_max,
                target_element_lower,
                target_element_upper,
                target_density_lower,
                target_density_upper,
                cell_size_lower, 
                cell_size_upper,
                bonding_probability_lower, 
                bonding_probability_upper,
                crossover_type,
                max_speed,
                squish_step_size,
                ensemble_type,
                ensemble_t_start,
                ensemble_t_end,
                ensemble_t_rate,
                run_steps,
                k,
                molecules,
                num_molecule,
                composition_num,
                reactions_num,
                rxns,
                ref_energy,
                weight_energy,
                weight_density,
                weight_element,
                exp_pdf,
                weight_pdf
            )


def get_element_percentage(skl_file, target_element):
    """
    Obtain the atomic percentage of a given target element by calling the external
    'element_pct.py' script as a subprocess.
    """
    element_percentage = subprocess.run(["element_pct.py", skl_file, target_element],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, check=True)
    return float(element_percentage.stdout.strip())


def read_last_value(filename):
    """
    Read the last non-empty line from a file and return it as a float.
    """
    with open(filename, "r") as f_values:
        lines = [line.strip() for line in f_values if line.strip()]
        return float (lines[-1])


def run_lammps_simulations(population, population_size, num_cores, target_element,
                           molecules, num_molecule, composition_num, reactions_num,
                           rxns, max_speed, squish_step_size, ensemble_type,
                           ensemble_t_start, ensemble_t_end, ensemble_t_rate,
                           run_steps, gen_count, exp_pdf=None):
    """
    This function will create the necessary input files for a lammps condensation
    and run lammps for each member of the Genetic Algorithm (GA)
    """
    # Creating seperate directories for generation and each member, each contains 
    # the lamps.in file. Lammps for each member is ran from the corresponding directory
    os.makedirs(f"generation_{gen_count}", exist_ok=True)
    os.chdir(f"generation_{gen_count}")
    for member in range(1, population_size + 1):
        os.makedirs(str(member), exist_ok=True)
        os.chdir(str(member))
        with open("condense.in", "w") as lamps_input_file:
            lamps_input_file.write(f"composition {composition_num}\n")
            for mol, num in zip(molecules, num_molecule):
                lamps_input_file.write(f"{mol} family1 {num}\n")
             
            lamps_input_file.write(f"\ncell_size {population[member - 1][2]}\n\n")
            lamps_input_file.write(f"max_speed {max_speed}\n")
            #lamps_input_file.write(f"condensation_rate {population[member - 1][3]}\n\n")
            lamps_input_file.write(f"stage {population[member - 1][3]} {squish_step_size} "
                                   f"{ensemble_type} {ensemble_t_start} {ensemble_t_end} "
                                   f"{ensemble_t_rate} {run_steps}\n\n")
            lamps_input_file.write(f"target_density {population[member - 1][1]}\n\n")
            
            lamps_input_file.write(f"reactions {reactions_num}\n")

            for rxn_index in range(reactions_num):
                lamps_input_file.write(f"{rxns[rxn_index][0]} {rxns[rxn_index][1]} "
                                       f"{rxns[rxn_index][2]} {rxns[rxn_index][3]} "
                                       f"{population[member - 1][4]}\n")

        # Wrap each individual in a try/except so that one failed member (e.g.
        # condense dying on a missing angle type, or LAMMPS crashing on a
        # neighbor list overflow) does not kill the entire EA run. Instead, the
        # member is assigned worst-case fitness values (1e99) so the EA treats
        # it as a bad individual and discards it naturally via selection/kill_rate.
        in_lammps = False

        # Write a header to stderr before running so that any crash output
        # that follows in deadmd.e is immediately identifiable by generation
        # and member number without having to count through the file.
        print(f"[DEAD-MD] Gen {gen_count} Member {member}: starting",
              file=sys.stderr, flush=True)

        try:
            subprocess.run("condense", shell=True, check=True)
            os.chdir("lammps")
            in_lammps = True
            #subprocess.run(f"srun -N 1 -c {num_cores} lmp < lammps.in", shell=True, check=True)
            subprocess.run(f"lmp < lammps.in", shell=True, check=True)
            #subprocess.run("mpirun -np 2 lmp < lammps.in", shell=True, check=True)

            print(f" I am in directory {member}\n")
            subprocess.run(
                "dump2skl -d dump.coarse -a lammps.dat -f=-1",
                shell=True, check=True)
            # Calculate elemental percentage.
            pct = get_element_percentage("olcao.skl", target_element)
            with open("element_evolve", "a") as f_elmnt:
                f_elmnt.write(f"{pct}\n")
            # Compute RPDF when experimental PDF has been loaded for
            # fitness evaluation. rpdf reads olcao.skl from the current
            # directory and writes rpdf.plot here (inside lammps/).
            # fitness_function reads rpdf.plot later to compute the
            # scale-optimised R-factor against the experimental G(r).
            if exp_pdf is not None:
                subprocess.run("rpdf", shell=True, check=True)

        except subprocess.CalledProcessError as e:
            print(f"[DEAD-MD] Gen {gen_count} Member {member}: CRASHED"
                  f" (command: {e.cmd}), assigning worst fitness.",
                  file=sys.stderr, flush=True)
            print(f"Member {member} failed (command: {e.cmd}),"
                  f" assigning worst fitness.")
            # fitness_function reads totE_evolve, density_evolve, and element_evolve
            # from each member's lammps/ directory. Write sentinel values so the EA
            # can continue rather than crashing on missing files.
            if in_lammps:
                # Failure happened inside lammps/ -- write sentinels here directly.
                for fname in ["totE_evolve", "density_evolve", "element_evolve"]:
                    with open(fname, "a") as f:
                        f.write("1e99\n")
            else:
                # Failure happened before entering lammps/ (e.g. condense crashed).
                os.makedirs("lammps", exist_ok=True)
                for fname in ["lammps/totE_evolve", "lammps/density_evolve",
                              "lammps/element_evolve"]:
                    with open(fname, "a") as f:
                        f.write("1e99\n")
        finally:
            # Always return to the top-level run directory regardless of
            # where the failure occurred.
            if in_lammps:
                os.chdir("../..")
            else:
                os.chdir("..")
    os.chdir("..")


def fitness_function(population, population_size, gen_count, ref_energy,
                     weight_energy=1.0, weight_density=1.0,
                     weight_element=1.0, exp_pdf=None, weight_pdf=0.0):
    """
    Compute a generation-independent fitness score for each individual.

    Three or four terms are summed (lower total is better), depending on
    whether an experimental PDF has been supplied:

      1. Energy term -- signed deviation from the user-supplied reference
         energy (ref_energy in deadmd.in), normalized by its magnitude
         so large absolute energies do not dominate the sum:
             energy_term = (E - ref_energy) / |ref_energy|
         The numerator is NOT wrapped in abs because energy is strictly
         directional: lower is always better, with no target to hit.
         A member worse than the reference (E > ref_energy) receives a
         positive penalty; one that beats the reference (E < ref_energy)
         receives a negative reward that lowers its total fitness.
         Taking abs of the numerator would incorrectly penalize a very
         stable structure as much as an equally-deviant unstable one.

      2. Density term -- relative absolute deviation of the simulated
         density from each member's own target density (population[i][1]),
         which was the value given to LAMMPS as input:
             density_term = |actual - target| / target
         abs IS used here because density has a specific target value:
         overshooting and undershooting are equally wrong. Using the
         member's own target keeps this term generation-independent.

      3. Element percentage term -- same form as the density term, using
         the member's own target element percentage (population[i][0]):
             element_term = |actual - target| / target
         Same reasoning: there is a specific desired percentage, so
         both directions of deviation are penalized symmetrically.

      4. PDF R-factor term (optional) -- scale-optimised R-factor between
         the simulated G(r) from rpdf and the experimental G(r) truncated
         to 10 Angstrom. The optimal scale factor s is found analytically:
             s        = (G_sim . G_exp) / |G_sim|^2
             pdf_term = sqrt( sum((s*G_sim - G_exp)^2) / sum(G_exp^2) )
         Using an optimal scale makes the comparison insensitive to the
         overall amplitude of rpdf's unnormalised output; only peak
         positions and relative heights drive the score. Active only when
         exp_pdf is not None and weight_pdf > 0.

    Crashed members (sentinel value 1e99) bypass all computation and
    receive a sentinel fitness so selection discards them naturally.

    The three terms are combined as a weighted sum:
        fitness = weight_energy  * energy_term
                + weight_density * density_term
                + weight_element * element_term
    All weights default to 1.0 (equal emphasis). Raising a weight makes
    the GA optimise that quantity more aggressively; setting a weight to
    0.0 removes that term from the selection pressure entirely.

    Args:
        population:      numpy array of genomes, shape (population_size, 5)
        population_size: number of individuals in the current generation
        gen_count:       current generation index
        ref_energy:      user-supplied reference total energy (float, read
                         from deadmd.in keyword ref_energy).  Should be
                         the energy of the best-known or target structure.
        weight_energy:   multiplicative weight for the energy term
                         (default 1.0; keyword weight_energy in deadmd.in).
        weight_density:  multiplicative weight for the density term
                         (default 1.0; keyword weight_density in deadmd.in).
        weight_element:  multiplicative weight for the element-percentage term
                         (default 1.0; keyword weight_element in deadmd.in).
        exp_pdf:         numpy array of experimental G(r), shape (1000,),
                         covering 0.01-10.00 Angstrom at 0.01 Angstrom
                         spacing. None when the PDF term is not requested.
        weight_pdf:      multiplicative weight for the PDF R-factor term
                         (default 0.0; keyword weight_pdf in deadmd.in).
                         The PDF term is active only when exp_pdf is not
                         None and weight_pdf > 0.

    Returns:
        fitnesses: list of floats, one per member (lower is better).
        As a side effect, writes generation_{gen_count}/fitness_scores_gen.dat
        containing one row per member with columns: member, rank, fitness.
        Rank 1 is the best (lowest fitness) and rank N is the worst -- the
        opposite of the rank-based selection method, where the best member
        receives the highest rank (N) to give it the greatest selection
        probability. The rank here is purely a human-readable label for
        inspecting generation output and has no effect on the algorithm.
    """
    SENTINEL = 1e99
    fitnesses = []
    use_pdf   = (exp_pdf is not None and weight_pdf > 0.0)

    for i in range(population_size):
        member_dir = str(i + 1)

        totE_file = (f"generation_{gen_count}/{member_dir}"
                     f"/lammps/totE_evolve")
        density_file = (f"generation_{gen_count}/{member_dir}"
                        f"/lammps/density_evolve")
        element_pct_file = (f"generation_{gen_count}/{member_dir}"
                            f"/lammps/element_evolve")

        total_energy    = read_last_value(totE_file)
        actual_density  = read_last_value(density_file)
        actual_elem_pct = read_last_value(element_pct_file)

        # Pass crashed members straight through without computing terms.
        if total_energy >= SENTINEL:
            fitnesses.append(SENTINEL)
            with open(f"generation_{gen_count}/{member_dir}"
                      f"/lammps/fitness_scores_gen", "a") as fitness_file:
                fitness_file.write(f"{SENTINEL}\n")
            continue

        # Each member carries its own target density and element percentage
        # in its genome. Using these as per-member references keeps the
        # density and element terms generation-independent.
        target_density  = population[i][1]
        target_elem_pct = population[i][0]

        energy_term  = ((total_energy - ref_energy)
                        / (abs(ref_energy) + 1e-12))
        density_term = (abs(actual_density - target_density)
                        / (target_density + 1e-12))
        element_term = (abs(actual_elem_pct - target_elem_pct)
                        / (target_elem_pct + 1e-12))

        # --- Optional PDF R-factor term ----------------------------------
        # Read rpdf.plot written by rpdf (inside lammps/) and compute a
        # scale-optimised R-factor against the experimental G(r). The
        # optimal scalar s is determined analytically:
        #   s = (G_sim . G_exp) / |G_sim|^2
        # so the comparison is insensitive to rpdf's unnormalised output
        # amplitude -- only peak positions and relative heights matter.
        # If rpdf.plot is missing or degenerate the member is failed with
        # SENTINEL so the EA discards it via selection/kill_rate.
        pdf_term = 0.0
        if use_pdf:
            rpdf_file = (f"generation_{gen_count}/{member_dir}"
                         f"/lammps/rpdf.plot")
            try:
                sim_data = np.loadtxt(rpdf_file)
                G_sim    = sim_data[:, 1]
                # A sum-of-squares of zero means rpdf produced no pairs,
                # which indicates a degenerate or crashed simulation.
                if np.dot(G_sim, G_sim) < 1e-12:
                    raise ValueError("all-zero G(r) in rpdf.plot")
                s        = (np.dot(G_sim, exp_pdf)
                            / np.dot(G_sim, G_sim))
                pdf_term = np.sqrt(
                    np.sum((s * G_sim - exp_pdf) ** 2)
                    / (np.sum(exp_pdf ** 2) + 1e-12)
                )
            except (FileNotFoundError, ValueError, OSError) as pdf_err:
                print(f"[DEAD-MD] PDF fitness error gen {gen_count}"
                      f" member {i + 1}: {pdf_err}",
                      file=sys.stderr, flush=True)
                fitnesses.append(SENTINEL)
                with open(f"generation_{gen_count}/{member_dir}"
                          f"/lammps/fitness_scores_gen", "a") as ff:
                    ff.write(f"{SENTINEL}\n")
                with open(f"generation_{gen_count}/{member_dir}"
                          f"/lammps/pdf_rfactor_evolve", "a") as ff:
                    ff.write(f"{SENTINEL}\n")
                continue
            # Write the R-factor so record_* functions can read it
            # without re-loading and re-computing.
            with open(f"generation_{gen_count}/{member_dir}"
                      f"/lammps/pdf_rfactor_evolve", "a") as ff:
                ff.write(f"{pdf_term}\n")

        # Normalize weights so their sum always equals the number of
        # active terms. At equal weights the factor is 1, reproducing
        # the original three-term behaviour exactly. When the PDF term
        # is active num_terms becomes 4 and the total_weight grows
        # accordingly, keeping the fitness magnitude stable.
        num_terms    = 4 if use_pdf else 3
        total_weight = (weight_energy + weight_density + weight_element
                        + (weight_pdf if use_pdf else 0.0))
        norm         = num_terms / (total_weight + 1e-12)

        fitness = norm * (weight_energy  * energy_term
                          + weight_density * density_term
                          + weight_element * element_term
                          + (weight_pdf * pdf_term if use_pdf else 0.0))
        fitnesses.append(fitness)

        with open(f"generation_{gen_count}/{member_dir}"
                  f"/lammps/fitness_scores_gen", "a") as fitness_file:
            fitness_file.write(f"{fitness}\n")

    # Write a single generation-level summary so all member scores and
    # their ranks are visible in one place without opening individual
    # member directories. Rank 1 is the best (lowest fitness). Crashed
    # members carry the sentinel value 1e99 and naturally rank last.
    sorted_indices = sorted(range(population_size), key=lambda i: fitnesses[i])
    ranks = [0] * population_size
    for rank, idx in enumerate(sorted_indices, start=1):
        ranks[idx] = rank

    fitness_scores_gen = f"generation_{gen_count}/fitness_scores_gen.dat"
    with open(fitness_scores_gen, "w") as summary_file:
        summary_file.write(f"{'member':<10} {'rank':<8} {'fitness':<20}\n")
        for i, score in enumerate(fitnesses):
            summary_file.write(
                f"{i + 1:<10} {ranks[i]:<8} {score:<20.6f}\n")

    return fitnesses


def record_best_member(gen_count, fitness_scores, exp_pdf=None):
    """
    Find the best-performing member in the current generation and append
    its generation number, member index, fitness, total energy, density,
    element percentage, and (when PDF is active) PDF R-factor as one line
    to best_by_fitness_profile_per_gen.dat in the top-level run directory.

    Args:
        gen_count:      current generation number
        fitness_scores: list of fitness values for all members
        exp_pdf:        experimental G(r) array or None; used only to
                        determine whether the pdf_rfactor column should
                        be written.
    """
    best_idx    = int(np.argmin(fitness_scores))
    best_member = best_idx + 1
    gen_dir     = f"generation_{gen_count}/{best_member}/lammps"

    total_energy = read_last_value(f"{gen_dir}/totE_evolve")
    density      = read_last_value(f"{gen_dir}/density_evolve")
    element_pct  = read_last_value(f"{gen_dir}/element_evolve")
    fitness      = fitness_scores[best_idx]

    pdf_rfactor = None
    if exp_pdf is not None:
        try:
            pdf_rfactor = read_last_value(
                f"{gen_dir}/pdf_rfactor_evolve")
        except (FileNotFoundError, ValueError):
            pdf_rfactor = float("nan")

    with open("best_by_fitness_profile_per_gen.dat", "a") as out_file:
        if gen_count == 1:
            header = (f"{'gen':<6} {'best_member':<13} {'fitness':<14}"
                      f" {'total_energy':<14} {'density':<10}"
                      f" element_pct")
            if exp_pdf is not None:
                header += "  pdf_rfactor"
            out_file.write(header + "\n")
        line = (f"{gen_count:<6} {best_member:<13} {fitness:<14.6f}"
                f" {total_energy:<14.6f} {density:<10.6f}"
                f" {element_pct:.6f}")
        if exp_pdf is not None:
            line += f"  {pdf_rfactor:.6f}"
        out_file.write(line + "\n")


def record_generation_averages(population_size, gen_count, fitness_scores,
                               exp_pdf=None):
    """
    Compute the average fitness, total energy, density, element percentage,
    and (when PDF is active) PDF R-factor across all members of the current
    generation and append one line to averages_per_gen.dat.

    Crashed members are assigned sentinel values of 1e99 by the simulation
    runner and are excluded from all averages so that a single crash does
    not corrupt the generation statistics. If every member crashed, 'nan'
    is written so matplotlib renders a visible gap rather than a spike.

    Args:
        population_size: number of members in the population
        gen_count:       current generation number
        fitness_scores:  list of fitness values for all members
        exp_pdf:         experimental G(r) array or None; controls whether
                         the pdf_rfactor column is collected and written.
    """
    # Any value at or above this threshold is a crash sentinel.
    SENTINEL = 1e90

    total_energies = []
    densities      = []
    element_pcts   = []
    pdf_rfactors   = []

    for member in range(1, population_size + 1):
        gen_dir = f"generation_{gen_count}/{member}/lammps"
        total_energies.append(
            read_last_value(f"{gen_dir}/totE_evolve"))
        densities.append(
            read_last_value(f"{gen_dir}/density_evolve"))
        element_pcts.append(
            read_last_value(f"{gen_dir}/element_evolve"))
        if exp_pdf is not None:
            try:
                pdf_rfactors.append(
                    read_last_value(f"{gen_dir}/pdf_rfactor_evolve"))
            except (FileNotFoundError, ValueError):
                pdf_rfactors.append(float("nan"))

    # Filter crashed members out of each quantity independently so a
    # partial crash still produces a meaningful average from survivors.
    valid_fitness   = [v for v in fitness_scores  if v < SENTINEL]
    valid_energies  = [v for v in total_energies  if v < SENTINEL]
    valid_densities = [v for v in densities       if v < SENTINEL]
    valid_elements  = [v for v in element_pcts    if v < SENTINEL]
    valid_rfactors  = ([v for v in pdf_rfactors
                        if v < SENTINEL and not np.isnan(v)]
                       if exp_pdf is not None else [])

    avg_fitness      = (sum(valid_fitness)   / len(valid_fitness)
                        if valid_fitness   else float("nan"))
    avg_total_energy = (sum(valid_energies)  / len(valid_energies)
                        if valid_energies  else float("nan"))
    avg_density      = (sum(valid_densities) / len(valid_densities)
                        if valid_densities else float("nan"))
    avg_element_pct  = (sum(valid_elements)  / len(valid_elements)
                        if valid_elements  else float("nan"))
    avg_pdf_rfactor  = (sum(valid_rfactors)  / len(valid_rfactors)
                        if valid_rfactors  else float("nan"))

    with open("averages_per_gen.dat", "a") as out_file:
        if gen_count == 1:
            header = (f"{'gen':<6} {'avg_fitness':<14}"
                      f" {'avg_total_energy':<18} {'avg_density':<13}"
                      f" avg_element_pct")
            if exp_pdf is not None:
                header += "  avg_pdf_rfactor"
            out_file.write(header + "\n")
        line = (f"{gen_count:<6} {avg_fitness:<14.6f}"
                f" {avg_total_energy:<18.6f} {avg_density:<13.6f}"
                f" {avg_element_pct:.6f}")
        if exp_pdf is not None:
            line += f"  {avg_pdf_rfactor:.6f}"
        out_file.write(line + "\n")


def record_best_quantities(population, population_size, gen_count,
                           fitness_scores, exp_pdf=None):
    """
    For each generation, find the member that is best by four independent
    criteria and append one line per criterion to its own dat file in the
    top-level run directory. The four criteria are:

      1. best_by_fitness_per_gen.dat -- lowest overall fitness score
      2. best_by_energy_per_gen.dat  -- lowest total energy (LJ potential)
      3. best_by_density_per_gen.dat -- actual density closest to the
                                        member's own genome target density
      4. best_by_element_per_gen.dat -- actual element percentage closest
                                        to the member's own genome target
                                        element percentage

    Each criterion may identify a different member. Recording them
    separately supports the planned fitness function redesign
    (DESIGN.md section 5) where density and element percentage are
    penalised by deviation from target rather than by raw value.

    Args:
        population:      2D numpy array of genome values
                           col 0: target element percentage
                           col 1: target density
                           col 2: cell size
                           col 3: condensation rate
                           col 4: bonding probability
        population_size: number of members in the population
        gen_count:       current generation number
        fitness_scores:  list of fitness values (one per member)
    """
    # Any value at or above this threshold is a crash sentinel, not a
    # real simulation result. If all members crashed, nan is written so
    # the plot shows a gap rather than a sentinel spike.
    SENTINEL = 1e90
    NAN      = float("nan")

    # Read actual simulation results for every member.
    total_energies = []
    densities      = []
    element_pcts   = []
    pdf_rfactors   = []

    for member in range(1, population_size + 1):
        gen_dir = f"generation_{gen_count}/{member}/lammps"
        total_energies.append(
            read_last_value(f"{gen_dir}/totE_evolve"))
        densities.append(
            read_last_value(f"{gen_dir}/density_evolve"))
        element_pcts.append(
            read_last_value(f"{gen_dir}/element_evolve"))
        if exp_pdf is not None:
            try:
                pdf_rfactors.append(
                    read_last_value(f"{gen_dir}/pdf_rfactor_evolve"))
            except (FileNotFoundError, ValueError):
                pdf_rfactors.append(SENTINEL)

    # --- Criterion 1: lowest fitness ---
    best_fit_idx    = int(np.argmin(fitness_scores))
    best_fit_member = best_fit_idx + 1
    all_crashed     = fitness_scores[best_fit_idx] >= SENTINEL

    with open("best_by_fitness_per_gen.dat", "a") as out_file:
        if gen_count == 1:
            out_file.write(
                f"{'gen':<6} {'member':<8} {'fitness':<14}"
                f" {'total_energy':<14} {'density':<12} element_pct\n"
            )
        if all_crashed:
            out_file.write(
                f"{gen_count:<6} {'nan':<8} {'nan':<14}"
                f" {'nan':<14} {'nan':<12} nan\n"
            )
        else:
            out_file.write(
                f"{gen_count:<6} {best_fit_member:<8}"
                f" {fitness_scores[best_fit_idx]:<14.6f}"
                f" {total_energies[best_fit_idx]:<14.6f}"
                f" {densities[best_fit_idx]:<12.6f}"
                f" {element_pcts[best_fit_idx]:.6f}\n"
            )

    # --- Criterion 2: lowest total energy ---
    best_energy_idx    = int(np.argmin(total_energies))
    best_energy_member = best_energy_idx + 1

    with open("best_by_energy_per_gen.dat", "a") as out_file:
        if gen_count == 1:
            out_file.write(
                f"{'gen':<6} {'member':<8} total_energy\n"
            )
        if total_energies[best_energy_idx] >= SENTINEL:
            out_file.write(f"{gen_count:<6} {'nan':<8} nan\n")
        else:
            out_file.write(
                f"{gen_count:<6} {best_energy_member:<8}"
                f" {total_energies[best_energy_idx]:.6f}\n"
            )

    # --- Criterion 3: actual density closest to genome target density ---
    density_deviations = [
        abs(densities[i] - population[i][1])
        for i in range(population_size)
    ]
    best_density_idx    = int(np.argmin(density_deviations))
    best_density_member = best_density_idx + 1

    with open("best_by_density_per_gen.dat", "a") as out_file:
        if gen_count == 1:
            out_file.write(
                f"{'gen':<6} {'member':<8} {'actual_density':<16}"
                f" {'target_density':<16} deviation\n"
            )
        if densities[best_density_idx] >= SENTINEL:
            out_file.write(
                f"{gen_count:<6} {'nan':<8} {'nan':<16}"
                f" {'nan':<16} nan\n"
            )
        else:
            out_file.write(
                f"{gen_count:<6} {best_density_member:<8}"
                f" {densities[best_density_idx]:<16.6f}"
                f" {population[best_density_idx][1]:<16.6f}"
                f" {density_deviations[best_density_idx]:.6f}\n"
            )

    # --- Criterion 4: actual element pct closest to genome target ---
    element_deviations = [
        abs(element_pcts[i] - population[i][0])
        for i in range(population_size)
    ]
    best_element_idx    = int(np.argmin(element_deviations))
    best_element_member = best_element_idx + 1

    with open("best_by_element_per_gen.dat", "a") as out_file:
        if gen_count == 1:
            out_file.write(
                f"{'gen':<6} {'member':<8} {'actual_pct':<12}"
                f" {'target_pct':<12} deviation\n"
            )
        if element_pcts[best_element_idx] >= SENTINEL:
            out_file.write(
                f"{gen_count:<6} {'nan':<8} {'nan':<12}"
                f" {'nan':<12} nan\n"
            )
        else:
            out_file.write(
                f"{gen_count:<6} {best_element_member:<8}"
                f" {element_pcts[best_element_idx]:<12.6f}"
                f" {population[best_element_idx][0]:<12.6f}"
                f" {element_deviations[best_element_idx]:.6f}\n"
            )

    # --- Criterion 5: lowest PDF R-factor (best structural match) ------
    # Active only when experimental PDF data was provided. Lower R-factor
    # means the simulated G(r) most closely matches the experimental one.
    if exp_pdf is not None:

        # Rank all members by pdf_rfactor (rank 1 = lowest = best).
        # Crashed/missing members carry the sentinel and rank last.
        # The ranking is computed first so the winner's rank is available
        # when writing both the per-generation and simulation-wide files.
        sorted_rf_idx = sorted(
            range(population_size),
            key=lambda i: pdf_rfactors[i]
        )
        rf_ranks = [0] * population_size
        for rank, idx in enumerate(sorted_rf_idx, start=1):
            rf_ranks[idx] = rank

        valid_rf = [(v, i) for i, v in enumerate(pdf_rfactors)
                    if v < SENTINEL and not np.isnan(v)]
        if valid_rf:
            best_pdf_val, best_pdf_idx = min(valid_rf, key=lambda x: x[0])
            best_pdf_member = best_pdf_idx + 1
            best_pdf_rank   = rf_ranks[best_pdf_idx]
        else:
            best_pdf_val    = None
            best_pdf_member = None
            best_pdf_rank   = None

        # Simulation-wide file: one row per generation, best member only.
        # After appending this generation's result, re-read all rows and
        # rewrite the file with cross-generation ranks so rank 1 always
        # identifies the generation with the lowest pdf_rfactor overall.
        sim_pdf_file = "best_by_pdf_per_gen.dat"
        new_row = (
            gen_count,
            best_pdf_member if best_pdf_member is not None else float("nan"),
            best_pdf_val    if best_pdf_val    is not None else float("nan"),
        )

        # Load all previously recorded rows (skip header line).
        prior_rows = []
        try:
            with open(sim_pdf_file, "r") as rf:
                lines = rf.readlines()
            for line in lines[1:]:          # lines[0] is the header
                parts = line.split()
                prior_rows.append((
                    int(parts[0]),
                    parts[1],               # member (may be "nan")
                    float(parts[3]),        # pdf_rfactor (col 3, skip rank)
                ))
        except FileNotFoundError:
            prior_rows = []

        all_rows = prior_rows + [new_row]

        # Rank generations by pdf_rfactor; nan/crashed entries rank last.
        def _sort_key(row):
            v = row[2] if not isinstance(row[2], float) else row[2]
            return (np.isnan(v), v)

        sorted_gens = sorted(range(len(all_rows)), key=lambda i: _sort_key(all_rows[i]))
        gen_ranks   = [0] * len(all_rows)
        for rank, idx in enumerate(sorted_gens, start=1):
            gen_ranks[idx] = rank

        with open(sim_pdf_file, "w") as out_file:
            out_file.write(
                f"{'gen':<6} {'member':<8} {'rank':<6} pdf_rfactor\n"
            )
            for idx, row in enumerate(all_rows):
                g, m, v = row
                if np.isnan(v):
                    out_file.write(
                        f"{g:<6} {'nan':<8} {'nan':<6} nan\n"
                    )
                else:
                    out_file.write(
                        f"{g:<6} {str(m):<8} {gen_ranks[idx]:<6}"
                        f" {v:.6f}\n"
                    )

        # Per-generation file: all members ranked by pdf_rfactor,
        # mirroring fitness_scores_gen.dat. Rank 1 = best (lowest).
        pdf_rfactor_gen = f"generation_{gen_count}/pdf_rfactor_gen.dat"
        with open(pdf_rfactor_gen, "w") as rf_file:
            rf_file.write(f"{'member':<10} {'rank':<8} pdf_rfactor\n")
            for i, rfactor in enumerate(pdf_rfactors):
                rf_file.write(
                    f"{i + 1:<10} {rf_ranks[i]:<8} {rfactor:.6f}\n"
                )


def parents_selection(selection, fitness_scores, kill_rate, population_size, k):
    """
    Select the top performing individuals based on fitness (lower fitness is better
    and means a fitter memberand the chosen selection method.
    Returns indices of slected parents.
    """
    num_parents = int(population_size - (population_size * kill_rate / 100))

    # Make sure number of parents is at least 2 and an even numebr to enable crossover.
    if num_parents < 2:
        num_parents = 2
    if num_parents % 2 != 0:
        num_parents += 1

    # Selection option 1 is Truncation selection
    if selection == 1:
        # Pair index with fitness, sort by fitness ascending
        ranked = sorted(
                enumerate(fitness_scores),
                key=lambda x: x[1])
    

        # Keep the top performers according to the kill_rate

        parents_indices = [indx for indx, score in ranked[:num_parents]]

    # Selection option 2 is rank selection
    elif selection == 2:
        # Rank individuals in an ascending order
        ranked = sorted(
                range(population_size),
                key=lambda i: fitness_scores[i]
        )        

        # assign linear ranks (best gets higher rank)
        ranks = np.array([population_size -i for i in range(population_size)], dtype=float)

        # Normalize to probabilities
        probabilities = ranks / np.sum(ranks)

        # Sample parents positions (with replacement)
        chosen_positions = np.random.choice(
                population_size,
                size=num_parents,
                replace=True,
                p=probabilities
        )

        # Map back to actual population indices
        parents_indices = [ranked[pos] for pos in chosen_positions]

    # Selection option 3 is tournament selection
    elif selection == 3:
        parents_indices = []

        for _ in range(num_parents):

            # Pick k random individuals to compete
            competitors = random.sample(range(population_size), k)

            # Select the fittest of that small tournament (lowest fitness)
            best = min(competitors, key=lambda i: fitness_scores[i])

            # add to parents list
            parents_indices.append(best)

    else:
        raise ValueError("selection must be truncation, rank or tournament")

    return parents_indices


def crossover_method(parent1, parent2, crossover_type):
    """
    Perform a crossover between two parents.

    Args:
        parent1 (list or np.ndarray): Genome of parent 1
        parent2 (list or np.ndarray): Genome of parent 2
    Returns:
    child1, child2 (same type as parents)
    """
    genome_length = len(parent1)

    # One point crossover
    if crossover_type == 1:

        crossover_point = random.randint(1, genome_length -1)

        child1 = np.concatenate((parent1[:crossover_point], parent2[crossover_point:]))
        child2 = np.concatenate((parent2[:crossover_point], parent1[crossover_point:]))
    
    # Two points crossover
    elif crossover_type == 2:
        crossover_point1, crossover_point2 = sorted(random.sample(range(1, genome_length), 2))
        child1 = np.concatenate(
                (parent1[:crossover_point1], parent2[crossover_point1:crossover_point2],
                 parent1[crossover_point2:])
                ) 
        child2 = np.concatenate(
                (parent2[:crossover_point1], parent1[crossover_point1:crossover_point2],
                 parent2[crossover_point2:])
                ) 



    return child1, child2


def mutation(child, mutation_rate, target_element_lower, target_element_upper,
             target_density_lower, target_density_upper, cell_size_lower,
             cell_size_upper, bonding_probability_lower, bonding_probability_upper):

    """
    Mutate a child genome by randomly changing genes with a giving mutation 
    probability.

    Args:
        child (list or np.ndarray): Genome to mutate
        mutation_rate (float): Probability of mutating each gene

    Returns:
        np.ndarray: Mutated genome
    """

    mutated_child = child.copy()

    # gene 0: target element percentage
    if random.random() < mutation_rate:
        mutated_child[0] = random.uniform(target_element_lower, target_element_upper)

    # gene 1: target density
    if random.random() < mutation_rate:
        mutated_child[1] = random.uniform(target_density_lower, target_density_upper)

    # gene 2: cell_size
    if random.random() < mutation_rate:
        mutated_child[2] = random.uniform(cell_size_lower, cell_size_upper)

    # gene 3: condensation rate
    if random.random() < mutation_rate:
        mutated_child[3] = random.uniform(0, 1)

    # gene 4: bonding probability
    if random.random() < mutation_rate:
        mutated_child[4]= random.uniform(bonding_probability_lower, bonding_probability_upper)

    return mutated_child


def create_offsprings(population, parents_indices, mutation_rate,
                      target_element_lower, target_element_upper,
                      target_density_lower, target_density_upper,
                      cell_size_lower, cell_size_upper,
                      bonding_probability_lower, bonding_probability_upper,
                      crossover_type, num_offspring):
    """
    Generate exactly num_offspring children from selected parents.

    Two distinct parents are picked at random from parents_indices each
    iteration, crossed over, and mutated. The loop repeats until at least
    num_offspring children exist, then the list is sliced to exactly that
    count. This guarantees the caller always receives a full generation
    worth of children regardless of kill_rate and elitism_rate settings.

    Picking randomly (rather than sequentially) gives every fit parent an
    equal mating probability and avoids pairing bias.

    Args:
        population (list): Current population, where each element is a
            genome.
        parents_indices (list): Indices of selected parents in the
            population.
        num_offspring (int): Exact number of children to produce. Should
            be population_size - elitism_count so that elites + children
            always reconstitute a full-size generation.

    Returns:
        list: Exactly num_offspring child genomes produced via crossover
              and mutation.
    """
    children = []

    while len(children) < num_offspring:
        # Pick two distinct parents at random from the fit pool.
        idx1, idx2 = random.sample(range(len(parents_indices)), 2)
        p1 = population[parents_indices[idx1]]
        p2 = population[parents_indices[idx2]]

        child1, child2 = crossover_method(p1, p2, crossover_type)

        child1 = mutation(child1, mutation_rate, target_element_lower,
                          target_element_upper, target_density_lower,
                          target_density_upper, cell_size_lower,
                          cell_size_upper, bonding_probability_lower,
                          bonding_probability_upper)
        child2 = mutation(child2, mutation_rate, target_element_lower,
                          target_element_upper, target_density_lower,
                          target_density_upper, cell_size_lower,
                          cell_size_upper, bonding_probability_lower,
                          bonding_probability_upper)

        children.append(child1)
        children.append(child2)

    return children[:num_offspring]


def compute_elitism_count(population_size, elitism_rate):
    """
    Return the number of elite individuals to carry into the next generation.

    Centralises the elitism logic so that create_offsprings (which needs
    to know how many slots are left for children) and new_generation (which
    fills those slots) always agree on the count.

    Args:
        population_size (int): Total number of individuals in the population.
        elitism_rate (float): Percentage of population preserved as elites.

    Returns:
        int: Number of elites (>= 1 when elitism_rate > 0, else 0).
    """
    count = int(elitism_rate * population_size / 100)
    # Ensure at least one elite survives when elitism is enabled.
    if elitism_rate > 0 and count == 0:
        count = 1
    return min(count, population_size)


def new_generation(population, fitness_scores, children, elitism_count):
    """
    Create the next generation population.

    Args:
        population (list): Current population (list of genomes)
        fitness_scores (list): Fitness score per individual
        children (list): Offspring genomes
        elitism_count (int): Number of top individuals to carry over
            unchanged. Compute this with compute_elitism_count() before
            calling.

    Returns:
        list: New generation
    """
    population_size = len(population)

    ranked_indices = sorted(
        range(population_size),
        key=lambda i: fitness_scores[i]
    )

    elites = [population[i] for i in ranked_indices[:elitism_count]]

    # ---- Fill the rest with children ----
    new_gen = elites.copy()

    for child in children:
        if len(new_gen) < population_size:
            new_gen.append(child)
        else:
            break

    return new_gen


def main():

    print_banner()

    (
        population, 
        population_size, 
        mutation_rate, 
        kill_rate, 
        elitism_rate,
        target_element, 
        num_cores, 
        selection, 
        gen_max,
        target_element_lower,
        target_element_upper,
        target_density_lower,
        target_density_upper,
        cell_size_lower, 
        cell_size_upper,
        bonding_probability_lower, 
        bonding_probability_upper,
        crossover_type,
        max_speed,
        squish_step_size,
        ensemble_type,
        ensemble_t_start,
        ensemble_t_end,
        ensemble_t_rate,
        run_steps,
        k,
        molecules,
        num_molecule,
        composition_num,
        reactions_num,
        rxns,
        ref_energy,
        weight_energy,
        weight_density,
        weight_element,
        exp_pdf,
        weight_pdf
    ) = initialize_population()

    gen_count = 1

    # -----GA loop excution -----
    while gen_count <= gen_max:
        print(f"---- Current generation is: {gen_count}/{gen_max} -----", flush=True)
        print(f"---- Population # in current gen is: {len(population)} ---", flush=True)
        run_lammps_simulations(population, population_size, num_cores,
                               target_element, molecules, num_molecule,
                               composition_num, reactions_num, rxns,
                               max_speed, squish_step_size, ensemble_type,
                               ensemble_t_start, ensemble_t_end,
                               ensemble_t_rate, run_steps, gen_count,
                               exp_pdf)
        fitness_scores = fitness_function(population, population_size,
                                          gen_count, ref_energy,
                                          weight_energy, weight_density,
                                          weight_element, exp_pdf,
                                          weight_pdf)
        record_best_member(gen_count, fitness_scores, exp_pdf)
        record_generation_averages(population_size, gen_count,
                                   fitness_scores, exp_pdf)
        record_best_quantities(population, population_size, gen_count,
                               fitness_scores, exp_pdf)
        print (fitness_scores)
        parents_indices = parents_selection(
            selection, fitness_scores, kill_rate, population_size, k)
        print(parents_indices)

        elitism_count = compute_elitism_count(population_size, elitism_rate)
        num_offspring = population_size - elitism_count

        children = create_offsprings(population, parents_indices,
                                     mutation_rate,
                                     target_element_lower,
                                     target_element_upper,
                                     target_density_lower,
                                     target_density_upper,
                                     cell_size_lower, cell_size_upper,
                                     bonding_probability_lower,
                                     bonding_probability_upper,
                                     crossover_type, num_offspring)
        # print(f"Children are {children}\n")
        population = new_generation(population, fitness_scores, children,
                                    elitism_count)
        # print(f"new population is {population}\n")

        gen_count += 1

        
if __name__ == "__main__":
        main()
