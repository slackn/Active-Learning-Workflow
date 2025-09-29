from ase.data import atomic_numbers 
from ase import Atoms
from ase.db import connect
from ase.ga.startgenerator import StartGenerator
from ase.ga.utilities import closest_distances_generator
from ase.ga.data import PrepareDB
import numpy as np
import argparse, yaml
from pathlib import Path

def create_db(cfg:dict,iteration):
    #Parameters
    n_atoms=cfg["initialization"]["n_atoms"]
    charge=cfg["initialization"]["charge"]
    n_to_generate=cfg["initialization"]["n_to_generate"]
    element=cfg["initialization"]["element"]
    db_name = f"data/iter{iteration:03d}/iter{iteration:03d}_{element}{n_atoms}_q{charge}.db"

    Z=atomic_numbers[element]

    # Compute bounding cube side length
    #side = 2 * (0.5 + ((3 * n_atoms) / (4 * np.pi * np.sqrt(2))) ** (1/3))
    side=5

    # 1. Define a fake slab (empty unit cell with no atoms)
    # This is the simulation box
    buffer=5    #extra vacuum space
    slab_side= side+buffer
    slab = Atoms(cell=[slab_side, slab_side, slab_side], pbc=False)

    # 2. Define what you're placing
    blocks = [(element, n_atoms)]

    # 3. Set minimum interatomic distances
    # Not sure if it works
    blmin = closest_distances_generator(atom_numbers=[Z], 
                                    ratio_of_covalent_radii=0.7)


    # Define the box in which to place atoms (inside the slab)
    #This s the placement box 
    origin = [(slab_side - side) / 2] * 3  # Center the box inside the slab
    box = [origin,
        [[side, 0, 0],
            [0, side, 0],
            [0, 0, side]]]

    db_path = Path(db_name)
    if db_path.exists():
        db_path.unlink()

    # 5. Prepare the database (overwrite if exists)
    d = PrepareDB(db_name, simulation_cell=slab, stoichiometry=[Z] * n_atoms)
        
    # 6. Create the StartGenerator
    sg = StartGenerator(
        slab=slab,
        blocks=blocks,
        blmin=blmin,
        box_to_place_in=box,
    )

    # 7. Generate and write the clusters
    for _ in range(n_to_generate):
        #print("Starting to create structure ", _)
        atoms = sg.get_new_candidate()
        #print("Start generator for structure ", _ , "is completed.")
        atoms.charge = charge
        atoms.info["charge"] = charge
        atoms.set_initial_charges([charge / n_atoms] * n_atoms)
        d.add_unrelaxed_candidate(atoms)
        print("Structure ", _ , " is sucessfully created.")

    return db_name

def main():
    ap= argparse.ArgumentParser(description="Create initial random population")
    ap.add_argument("--config", "-c", required=True, help="Path to project config.yaml")
    ap.add_argument("--iter", "-i", required=True, type=int, help="Iteration number to save into correct folder")
    args=ap.parse_args()
    cfg=yaml.safe_load(open(args.config))
    iteration=args.iter
    create_db(cfg, iteration)

if __name__ == "__main__":
    main()