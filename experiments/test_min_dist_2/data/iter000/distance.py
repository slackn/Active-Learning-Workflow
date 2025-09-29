import numpy as np
from ase.io import read

# --- Read all frames from the .extxyz file ---
frames = read("selected_for_dft.extxyz", index=":")  # ':' = all frames

for frame_idx, atoms in enumerate(frames):
    coords = atoms.get_positions()  # shape (N_atoms, 3)
    num_atoms = len(coords)
    
    # Compute distance matrix
    distances = np.zeros((num_atoms, num_atoms))
    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            dist = np.linalg.norm(coords[i] - coords[j])
            distances[i, j] = dist
            distances[j, i] = dist
    
    # Minimum interatomic distance (ignore zeros on diagonal)
    min_dist = np.min(distances[np.triu_indices(num_atoms, k=1)])
    
    # Print results
    np.set_printoptions(precision=4, suppress=True)
    print(f"\nFrame {frame_idx}:")
    print("Distance matrix (Å):")
    print(distances)
    print(f"Minimum interatomic distance: {min_dist:.4f} Å")
