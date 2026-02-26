import sys
from ase.calculators.emt import EMT
from ase.io import read

# 1. Structural relaxation. ##################################################

# Check if a POSCAR file is provided
if len(sys.argv) != 2:
    print("Usage: python get_emt.py <POSCAR_file>")
    sys.exit(1)

# Input file path
poscar_file = sys.argv[1]

try:
    # Setup calculator:
    ase_calculator = EMT()

    # Read the POSCAR file
    slab_final = read(poscar_file)
    slab_final.set_calculator(ase_calculator)

    # Calculate energy (eV) and forces (eV/Å)
    energy = slab_final.get_potential_energy()
    forces = slab_final.get_forces()

    # Write energy to a text file
    with open('energy_output.txt', 'w') as energy_file:
        energy_file.write(f"{energy}\n")

    # Write forces to a text file
    with open('forces_output.txt', 'w') as forces_file:
        for atom_idx, force in enumerate(forces):
            forces_file.write(f"{force[0]:.6f} {force[1]:.6f} {force[2]:.6f}\n")

    print("Energy and forces have been successfully written to output files.")

except Exception as e:
    print(f"Error: {e}")

