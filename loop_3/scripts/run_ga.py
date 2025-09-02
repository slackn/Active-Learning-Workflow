import yaml 
import argparse
from random import random
from scripts.committee_calc import CommitteeCalculator
from ase.optimize import BFGS
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.data import DataConnection
from ase.ga.offspring_creator import OperationSelector
from ase.ga.population import Population
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.standardmutations import MirrorMutation, RattleMutation, PermutationMutation, StrainMutation
from ase.ga.utilities import closest_distances_generator, get_all_atom_types
from ase.data import atomic_numbers
from ase.io import write
from pathlib import Path


def run_ga(cfg:dict, iteration):

    # Parameters
    n_atoms=cfg["initialization"]["n_atoms"]
    charge=cfg["initialization"]["charge"]
    n_to_generate=cfg["initialization"]["n_to_generate"]
    element=cfg["initialization"]["element"]
    db_name = f"data/iter{iteration:03d}/iter{iteration:03d}_{element}{n_atoms}_q{charge}.db"
    offsprings= cfg["ga"]["offsprings"]
    mutation_probability= cfg["ga"]["mutation_prob"]
    fmax=cfg["ga"]["fmax"]
    opt_steps= cfg["ga"]["opt_steps"]
    population_size=cfg["initialization"]["n_to_generate"]
    


    Z=atomic_numbers[element]


    # Connect to database
    db= DataConnection(db_name)
    atom_numbers= [Z] * n_atoms

    # Define geometric constraints
    slab= db.get_slab()
    all_atom_types=get_all_atom_types(slab, atom_numbers)
    blmin= closest_distances_generator(all_atom_types, ratio_of_covalent_radii=0.7)

    # Comparator: to get distinct canidadates
    comparator= InteratomicDistanceComparator(
        n_top=n_atoms,
        pair_cor_cum_diff=0.015,
        pair_cor_max=0.7,
        #dE=0.02,
        mic= False
    )

    # Operators
    pairing= CutAndSplicePairing(slab=slab, n_top=n_atoms, blmin=blmin)
    mutations= OperationSelector(
        [1.0, 1.0],
        [
            MirrorMutation(blmin,n_atoms),
            RattleMutation(blmin, n_atoms)
        ]
        #StrainMutation
    )

    # Set MACE model as calculator 

    committee_calc= CommitteeCalculator(
        iteration=iteration,
        use_forces=True
    )

    # Relax initial random structures

    while db.get_number_of_unrelaxed_candidates()>0:
        atoms= db.get_an_unrelaxed_candidate()
        confid=atoms.info.get('confid', 'N/A')

        atoms.set_calculator(committee_calc)    # this gives a future warning
        atoms.calc=committee_calc
        print(f"[GA] Relaxing initial random structure confid={confid}")
        dyn=BFGS(atoms,logfile="opt_initial.log")
        dyn.run(fmax, opt_steps)

        E=atoms.get_potential_energy()
        atoms.info['key_value_pairs']['raw_score'] = -E
        db.add_relaxed_step(atoms)
        
        energy= atoms.calc.results["energy"]
        print(f"[Relax] Energy from calculator = {energy}")

    # Build population from relaxed structures
    population= Population(
        data_connection=db,
        population_size=population_size,
        comparator=comparator
    )

    print(f"\n[Population] Current population size: {len(population.pop)}")


    # Main GA loop
    for i in range(offsprings):
        print(f"\n[GA] Starting candidate {i + 1}/{offsprings}")
        parent1, parent2 = population.get_two_candidates()
        # --- Print parent info ---
        confid1 = parent1.info.get('confid', 'N/A')
        confid2 = parent2.info.get('confid', 'N/A')
        print(f"[GA] Selected parents:")
        print(f"   Parent 1: confid={confid1}")
        print(f"   Parent 2: confid={confid2}")

        child, description = pairing.get_new_individual([parent1, parent2])
        if child is None:
            print("[GA] Pairing failed. Skipping.")
            continue
        db.add_unrelaxed_candidate(child, description)

        if random() < mutation_probability:
            mutated_child, desc = mutations.get_new_individual([child])
            if mutated_child is not None:
                db.add_unrelaxed_step(mutated_child, desc)
                child = mutated_child
                print("[GA] Mutation applied.")

        child.set_calculator(committee_calc)    # this gives a future warning
        child.calc=committee_calc
        dyn=BFGS(child, logfile="child_opt.log")
        print('Offspring relaxation starts.')
        dyn.run(fmax, opt_steps)

        E=child.get_potential_energy()
        child.info['key_value_pairs']['raw_score'] = -E
    
        db.add_relaxed_step(child)
        population.update()


    candidates_list=list(population.pop)


    # Sort relaxed candidates by decreasing sigma_E_pa
    candidates_list.sort(
        key=lambda a:float(a.info.get("key_value_pairs", {}).get("sigma_E_pa", float("nan"))
        ),
        reverse=True
    )
    # Print confid and energy uncertainty for each structure       
    for structure in candidates_list:
        confid= structure.info.get('confid', 'N/A')
        kv=structure.info.get("key_value_pairs",{})
        sigma_E_pa= float(kv.get("sigma_E_pa", float("nan")))
        sigma_F_mean=float(kv.get("sigma_F_mean", float("nan")))
        raw_score=kv.get("raw_score")
        print(f"confid={confid}, raw score={raw_score}, energy uncertainty={sigma_E_pa}, force uncertainty={sigma_F_mean}")


    n_dft=cfg["ga"]["n_dft"]
    # Save highest n_dft number of structures 
    n_dft = min(n_dft, len(candidates_list))    #candidates list might be smaller than n_dft
    selected = candidates_list[:n_dft]

    out_xyz=f"data/iter{iteration:03d}/selected_for_dft.extxyz"
    Path(out_xyz).parent.mkdir(parents=True, exist_ok=True)

    write(out_xyz, selected, format="extxyz")  # classic XYZ, no structured metadata
    print(f"Wrote {len(selected)} structures → {out_xyz}")


def main():
    ap= argparse.ArgumentParser(description="Suggest structure for dft labeling via genetic algorithm")
    ap.add_argument("--config", "-c", required=True, help="Path to project config.yaml")
    ap.add_argument("--iter","-i", required=True,type=int, help="Iteration index for genetic algorithm")
    args=ap.parse_args()
    cfg=yaml.safe_load(open(args.config))
    iteration=args.iter
    run_ga(cfg,iteration)
    
if __name__ == "__main__":
        main()