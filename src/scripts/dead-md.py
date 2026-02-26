#!/usr/bin/env python3


import os
import numpy as np
import random
import subprocess


def initialize_population():
    """ 
    Reads the input file dead-md.in, and then based on that it will create 
    an initial population of N LAMMPS simulations with their necessary 
    files (data files, input files, reaction templates).

    Args:
        None

    Returns:
        Population_size, elitism_rate, 
    """
    population = []
    molecules = []
    molecule_family = []
    num_molecule = []
    rxns = []

# Read the dead-md.in input file
    with open("dead-md.in", "r") as input_file:
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
                mutation_rate = int(words[1])
            elif (words[0].lower() == "kill_rate"):
                kill_rate = int(words[1])
            elif (words[0].lower() == "crossover"):
                crossover = int(words[1])
            elif (words[0].lower() == "selection"):
                selection_method = words[1]
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
        # The sequence of probabilites in the array is in respect to the seuqence
        #  of the reation you list in the dead-md.in or lamps.in
        population = np.empty((population_size, 5))
        for i in range(population_size):
            population[i][0] = random.uniform(target_element_lower, target_element_upper)
            population[i][1] = random.uniform(target_density_lower, target_density_upper)
            population[i][2] = random.uniform(cell_size_lower, cell_size_upper)
            population[i][3] = random.uniform(0.4, 0.99)
            population[i][4] = random.uniform(bonding_probability_lower, 
                    bonding_probability_upper)
        print (population)

    # Creating seperate directories for each member each contains the lamps.in
    #  file. Lammps for each member is ran from the corresponding directory
    for member in range(1, population_size + 1):
        os.makedirs(str(member), exist_ok=True)
        os.chdir(str(member))
        with open("condense.in", "w") as lamps_input_file:
            lamps_input_file.write(f"composition {composition_num}\n")
            for mol, num in zip(molecules, num_molecule):
                lamps_input_file.write(f"{mol} family1 {num}\n")
             
            lamps_input_file.write(f"\ncell_size {population[member - 1][2]}\n\n")
            lamps_input_file.write(f"condensation_rate {population[member - 1][3]}\n\n")
            lamps_input_file.write(f"target_density {population[member - 1][1]}\n\n")
            
            lamps_input_file.write(f"reactions {reactions_num}\n")

            for rxn_index in range(reactions_num):
                lamps_input_file.write(f"{rxns[rxn_index][0]} {rxns[rxn_index][1]} "
                                       f"{rxns[rxn_index][2]} {rxns[rxn_index][3]} "
                                       f"{population[member - 1][4]}\n")
            
        
        os.chdir("..")

    return (population, 
            population_size,
            mutation_rate,
            kill_rate,
            elitism_rate, 
            target_element,
            num_cores, crossover, selection_method, gen_max)



def get_element_percentage(skl_file, element):
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


def run_lammps_simulations(population_size, num_cores):
    """
    This function will create the necessary input files for a lammps condensation
    and run lammps for each member of the Genetic algorithm (GA)
    """
    for member in range(1, population_size + 1):
        os.chdir(str(member))
        subprocess.run("condense", shell=True, check=True)
        os.chdir("lammps")
        subprocess.run(f"srun -N 1 -c {num_cores} lmp < lammps.in", shell=True, check=True)
        #subprocess.run("mpirun -np 2 lmp < lammps.in", shell=True, check=True)

        print(f" I am in directory {member}\n")
        subprocess.run("dump2skl -d dump.coarse -a lammps.dat -f=-1", shell=True, check=True)
        # Calculate elemental percentage
        pct = get_element_percentage("olcao.skl", target_element)
        with open("element_evolve", "a") as f_elmnt:
            f_elmnt.write(f"{pct}\n")
        os.chdir("../..")


        



def fitness_function(population_size):
    """
    Compute fitness score for each individual.
    Fitness (for now) = total_energy + density + element_percentage

    Returns:
        fitness_scores: list of floats (lower is better)
    """
    fitnesses = []
    for i in range(population_size):
        member_dir = str(i + 1)

        totE_file = f"{member_dir}/lammps/totE_evolve"
        density_file = f"{member_dir}/lammps/density_evolve"
        element_pct_file = f"{member_dir}/lammps/element_evolve"

        total_energy = read_last_value(totE_file)
        density = read_last_value(density_file)
        element_pct  = read_last_value(element_pct_file)

        fitness = total_energy + density + element_pct
        fitnesses.append(fitness)
        # Storing the fitness scores.
        with open(f"{member_dir}/lammps/fitness_scores_gen", "a") as fitness_file:
            fitness_file.write(f"{fitness}\n")
    return fitnesses






def parents_selection(fitness_scores, kill_rate, population_size):
    """
    Select the top performing individuals based on fitness.
    Returns indices of slected parents.
    """
    num_parents = int(population_size * kill_rate / 100)

    # Make sure number of parents is at least 2 and even to enable crossover.
    if num_parents < 2:
        num_parents = 2
    if num_parents % 2 != 0:
        num_parents += 1

    # Pair index with fitness, sort by fitness descending

    ranked = sorted(enumerate(fitness_scores), key=lambda x: x[1], reverse=True)
    
    # Keep the top performers according to the kill_rate

    indices = [indx for indx, score in ranked[:num_parents]]

    return indices




def crossover_method(parent1, parent2):
    """
    Perform a crossover between between two parents.

    Args:
        parent1 (list or np.ndarray): Genome of parent 1
        parent2 (list or np.ndarray): Genome of parent 2
    Returns:
    child1, child2 (same type as parents)
    """
    genome_length = len(parent1)

    crossover_point = random.randint(1, genome_length -1)

    child1 = np.concatenate((parent1[:crossover_point], parent2[crossover_point:]))
    child2 = np.concatenate((parent2[:crossover_point], parent1[crossover_point:]))

    return child1, child2


def create_offsprings(population, parents_indices):
    """
    Generate offspring from selected parents using crossover.

    This function pairs parents two-by-two (based on their indices),
    applies the crossover operator, and produces a new set of child
    genomes. Parents themselves are not modified.

    Args:
        population (list): Current population, where each element is a genome.
        parent_indices (list): Indices of selected parents in the population.

    Returns:
        list: A list of child genomes produced via crossover.
              The number of children equals the number of parents
              (2 parents → 2 children).
    """
    children = []

    for i in range(0, len(parents_indices) - 1, 2):
        p1 = population[parents_indices[i]]
        p2 = population[parents_indices[i + 1]]

        child1, child2 = crossover_method(p1, p2)

        children.append(child1)
        children.append(child2)

    return children



def new_generation(population, fitness_scores, children, elitism_rate):
    """
    Create the next generation population.

    Args:
        population (list): Current population (list of genomes)
        fitness_scores (list): Fitness score per individual
        children (list): Offspring genomes
        elitism_rate (float): Fraction of population to keep as elites

    Returns:
        list: New generation
    """
    population_size = len(population)

    # ---- Elitism ----
    elitism_count = int(elitism_rate * population_size / 100)

    if elitism_rate > 0 and elitism_count == 0:
        elitism_count = 1

    elitism_count = min(elitism_count, population_size)

    ranked_indices = sorted(
        range(population_size),
        key=lambda i: fitness_scores[i],
        reverse=True
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


population, population_size, mutation_rate, kill_rate, elitism_rate, target_element, num_cores, crossover, slelction_method, gen_max = initialize_population()
gen_count = 1

# -----GA loop excution -----
while gen_count <= gen_max:
    run_lammps_simulations(population_size, num_cores)
    fitness_scores = fitness_function(population_size)
    print (fitness_scores)
    parents_indices = parents_selection(fitness_scores, kill_rate, population_size)
    print (parents_indices)
    children = create_offsprings(population, parents_indices)
    print(f"Children are {children}\n")
    population = new_generation(population, fitness_scores, children, elitism_rate)
    print(f"new population is {population}\n")

    gen_count += 1
