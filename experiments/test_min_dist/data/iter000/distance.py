import numpy as np

# Function to read xyz file and extract coordinates
def read_xyz(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    num_atoms = int(lines[0].strip())
    coords = []
    for line in lines[2:2 + num_atoms]:  # skip the first two lines
        parts = line.split()
        x, y, z = map(float, parts[1:4])
        coords.append([x, y, z])
    return np.array(coords)

# Function to compute distances between all pairs of atoms
def compute_distances(coords):
    num_atoms = len(coords)
    distances = np.zeros((num_atoms, num_atoms))
    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            dist = np.linalg.norm(coords[i] - coords[j])
            distances[i, j] = dist
            distances[j, i] = dist  # symmetric
    return distances

# Main
filename = "selected_for_dft.extxyz"
coords = read_xyz(filename)
distances = compute_distances(coords)

# Print distance matrix
np.set_printoptions(precision=4, suppress=True)
print("Distance matrix (Å):")
print(distances)
