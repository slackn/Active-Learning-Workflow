import shutil
from pathlib import Path
from ase.io import read, write
from ase.calculators.turbomole import Turbomole
from ase.units import Hartree, Bohr
import numpy as np

# ----- User settings -----
input_file = "new_dataset.xyz"  # your multi-frame .xyz/.extxyz
output_file = "dataset_relabelled.extxyz"
params = {
    'total charge': -1,
    'multiplicity': 1,
    'scf iterations': 1000,
    'basis set name': 'def-SV(P)',
    'density functional': 'pbe'

}

if shutil.which("define") is None:
    print("WARNING: Turbomole binaries not found in PATH. Input files will be generated but calculation won't run.")

# ----- Read all frames -----
frames = read(input_file, index=":")  # list of ASE Atoms objects

def get_dft_calculator():
    return Turbomole(**params)

# ----- Loop over frames, attach calculator, recompute -----
for i, atoms in enumerate(frames):
    print(f"Frame {i} is being relaxed. ")
    #calc = Turbomole(**tm_params)
    #atoms.set_calculator(calc)
    atoms.calc=get_dft_calculator()

    energy_eV = atoms.get_potential_energy()
    forces_eV_A = atoms.get_forces()

    energy_Ha = energy_eV / Hartree
    forces_Ha_Bohr = forces_eV_A / Hartree * Bohr

    # Save computed values as arrays in Atoms object (optional)
    atoms.info["energy_eV"] = energy_eV
    atoms.arrays["forces_eV_A"] = forces_eV_A
    atoms.arrays["forces_Ha_Bohr"] = forces_Ha_Bohr

    # Print summary for this frame
    print(f"\nFrame {i}: Energy = {energy_eV:.8f} eV ({energy_Ha:.10f} Ha)")
    print("# Forces (eV/Å) / (Ha/Bohr):")
    for sym, (fx, fy, fz), (Fx, Fy, Fz) in zip(
        [a.symbol for a in atoms], forces_eV_A, forces_Ha_Bohr
    ):
        print(f"{sym:>2s}  {fx: .6e} {fy: .6e} {fz: .6e}   ||  {Fx: .6e} {Fy: .6e} {Fz: .6e}")

# ----- Save all frames to a new extxyz -----
write(output_file, frames)
print(f"\nAll frames recomputed and saved to {output_file}")
