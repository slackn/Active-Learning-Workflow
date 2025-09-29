import re

input_file = "dataset_iter000.xyz"   # replace with your file
output_file = "output_mace.xyz"

with open(input_file, "r") as f:
    lines = f.readlines()

with open(output_file, "w") as f_out:
    i = 0
    while i < len(lines):
        # Skip empty lines
        if lines[i].strip() == "":
            i += 1
            continue

        # Number of atoms
        num_atoms = int(lines[i].strip())
        f_out.write(f"{num_atoms}\n")
        i += 1

        # Properties line
        props_line = lines[i].strip()
        energy_match = re.search(r"energy_eV=([-\d.]+)", props_line)
        charge_match = re.search(r"charge=([-\d.]+)", props_line)
        energy_eV = float(energy_match.group(1)) if energy_match else 0.0
        charge = int(float(charge_match.group(1))) if charge_match else 0

        # Write MACE header
        f_out.write(f"REF_energy={energy_eV:.16f} charge={charge} Properties=species:S:1:pos:R:3:REF_forces:R:3\n")
        i += 1

        # Write atom positions and forces, nicely aligned
        for _ in range(num_atoms):
            atom_line = lines[i].strip().split()
            symbol = atom_line[0]
            x, y, z = map(float, atom_line[1:4])
            first=map(float, atom_line[1])
            print(first)
            #print(x,y,z)
            #print(map(float, atom_line[1:4]))
            #fx, fy, fz = map(float, atom_line[7:9])  # forces_eV_A

            # Use fixed width formatting to match reference
            f_out.write(f"{symbol:<2}  "
                        f"{x:>22.16f}  {y:>22.16f}  {z:>22.16f}  "
                        f"{fx:>22.16f}  {fy:>22.16f}  {fz:>22.16f}\n")
            i += 1
