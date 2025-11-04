
import subprocess
import os
import random
from deap import base, creator, tools

# --- SETTINGS ---
N_POP = 10
N_GEN = 5
LAMMPS_EXEC = "lmp_serial"
CONDENSE_SCRIPT = "./condense"
DUMP_DIR = "logs"

os.makedirs(DUMP_DIR, exist_ok=True)

# --- EA SETUP ---
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()
toolbox.register("attr_param", random.random)
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_param, n=1)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

toolbox.register("mate", tools.cxBlend, alpha=0.5)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=0.2, indpb=0.5)
toolbox.register("select", tools.selTournament, tournsize=3)

def generate_structure(tag):
    print(f"Generating structure: {tag}")
    subprocess.run([CONDENSE_SCRIPT, tag])
    run_lammps(tag)

def run_lammps(tag):
    input_file = f"lammps_{tag}.in"
    log_file = os.path.join(DUMP_DIR, f"log_{tag}.lammps")
    print(f"Running LAMMPS for: {tag}")
    subprocess.run([LAMMPS_EXEC, "-in", input_file], stdout=open(log_file, "w"), stderr=subprocess.STDOUT)

def evaluate(tag):
    log_file = os.path.join(DUMP_DIR, f"log_{tag}.lammps")
    try:
        with open(log_file) as f:
            for line in f:
                if "Total Energy" in line:
                    return float(line.split()[-1]),
    except:
        return 1e6,
    return 1e6,

def main():
    pop = toolbox.population(n=N_POP)

    for i, ind in enumerate(pop):
        tag = f"gen0_indiv{i}"
        generate_structure(tag)
        ind.fitness.values = evaluate(tag)

    for gen in range(1, N_GEN + 1):
        print(f"\n--- Generation {gen} ---")
        offspring = toolbox.select(pop, len(pop))
        offspring = list(map(toolbox.clone, offspring))

        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < 0.7:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values

        for mutant in offspring:
            if random.random() < 0.3:
                toolbox.mutate(mutant)
                del mutant.fitness.values

        for i, ind in enumerate(offspring):
            if not ind.fitness.valid:
                tag = f"gen{gen}_indiv{i}"
                generate_structure(tag)
                ind.fitness.values = evaluate(tag)

        pop[:] = offspring

    best = tools.selBest(pop, 1)[0]
    print("\nBest individual:", best)

if __name__ == "__main__":
    main()
