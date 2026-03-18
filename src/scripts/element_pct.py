#!/usr/bin/env python3

import sys

def calculate_element_percentage(filename, target_element):
    """
    Calculate the atomic percentage of a given element in an OLCAO .skl structure file.

    The function parses the `fractional` section of an OLCAO .skl file, counts how many
    atoms belong to the specified element, and computes its percentage relative to the
    total number of atoms declared in the file.

    Parameters
    ----------
    filename : str
        Path to the OLCAO .skl file containing atomic fractional coordinates.
    target_element : str
        Chemical symbol of the element to evaluate (e.g. "h", "c", "si", "n").
        The comparison is case-insensitive.

    Returns
    -------
    float
        Atomic percentage of the target element in the structure.

    Raises
    ------
    ValueError
        If the total atom count cannot be determined from the file.
    FileNotFoundError
        If the input file does not exist.

    """
    reading_atoms = False
    element_count = 0
    total_atoms = 0

    with open(filename, 'r') as file:
        for line in file:
            stripped = line.strip().lower()

            if stripped.startswith("fractional"):
                parts = stripped.split()
                if len(parts) == 2 and parts[1].isdigit():
                    total_atoms = int(parts[1])
                    reading_atoms = True
                continue

            if reading_atoms and (stripped == "" or stripped.startswith("space")):
                break

            if reading_atoms:
                parts = stripped.split()
                if len(parts) >= 4:
                    atom_type = parts[0]
                    if atom_type == target_element.lower():
                        element_count += 1

    if total_atoms == 0:
        raise ValueError("Total atom count not found or is zero.")

    percentage = 100.0 * element_count / total_atoms
    return percentage

# Command line usage
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: get_element_percent.py <filename> <element>")
        sys.exit(1)

    filename = sys.argv[1]
    element = sys.argv[2]

    try:
        pct = calculate_element_percentage(filename, element)
        print(f"{pct:.2f}")
        #print(f"{element.upper()} content: {pct:.2f}%")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
