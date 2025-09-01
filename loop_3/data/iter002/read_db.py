""" from ase.ga.data import DataConnection

db = DataConnection("Na8_q0.db")  # adjust path/name

# Get relaxed candidates
for atoms in db.get_all_relaxed_candidates():
    confid = atoms.info.get("confid", "N/A")
    energy = atoms.get_potential_energy() if atoms.calc else "N/A"
    print(f"confid={confid}, energy={energy}")
 """

from ase.ga.data import DataConnection

db_name = "Na8_q0.db"  # adjust path
db = DataConnection(db_name)

# Get all relaxed candidates
relaxed_list = db.get_all_relaxed_candidates()
print(f"Number of relaxed candidates: {len(relaxed_list)}")

for atoms in relaxed_list:
    print("===")
    print("confid:", atoms.info.get("confid", "N/A"))
    print("generation:", atoms.info.get("generation", "N/A"))
    print("parent_ids:", atoms.info.get("parent_ids", "N/A"))
    print("raw_score:", atoms.info.get("key_value_pairs", {}).get("raw_score", "N/A"))
    print("sigma_E_pa:", atoms.info.get("key_value_pairs", {}).get("sigma_E_pa", "N/A"))
    print("sigma_F_mean:", atoms.info.get("key_value_pairs", {}).get("sigma_F_mean", "N/A"))
    print("other info keys:", list(atoms.info.keys()))
    print("positions:", atoms.get_positions())
    print("chemical symbols:", atoms.get_chemical_symbols())
    print("number of atoms:", len(atoms))
    if atoms.calc:
        try:
            print("energy:", atoms.get_potential_energy())
        except:
            print("energy: N/A")
    else:
        print("calculator: None")
