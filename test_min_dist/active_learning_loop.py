import argparse
import yaml 

from scripts.bootstrap import run_bootstrap
from scripts.train_mace import train_ensemble_for_iteration
from scripts.calc_mean_error import compute_mean_test_mae_for_iteration
from scripts.create_db import create_db
from scripts.run_ga import run_ga
from scripts.submit_dft import submit_dft
from scripts.merge import merge_datasets
from pathlib import Path
import matplotlib.pyplot as pyplot
from ase.io import read


# Parse arguments in command 
parser= argparse.ArgumentParser()
parser.add_argument("--config", "-c", required=True, help="Path to config file")
args= parser.parse_args()

# Load config file
with open(args.config) as f:
    cfg= yaml.safe_load(f)

# Use config
#print("Config loaded:", cfg)
#print("Learning rate: ",cfg["training"]["lr"])
#print("Batch size: ",cfg["training"]["batch_size"])

first_average_errors=[]
second_average_errors=[]


for it in range(0, cfg["active_learning"]["iterations"]):
    print(it)
    
    #Prepare bootstrapped training datasets
    manifest=run_bootstrap(cfg, iteration=it)

    manifest_path=Path(manifest["outdir"]) / "manifest.json"
    # Train ensemble for this iteration 
    train_ensemble_for_iteration(cfg, manifest_path)
    # Compute mean test MAE after training
    first_mae, second_mae = compute_mean_test_mae_for_iteration(it)
    first_average_errors.append(first_mae)
    second_average_errors.append(second_mae)
    print("Average mae per atom on test data in stage one: ", first_mae)
    print("Average mae per atom on test datain stage two: ", second_mae)
    # Create initial random population
    create_db(cfg, it)
    # Run genetic algorithm and select uncertain candidates
    run_ga(cfg, it)
    # Submit dft labeling 
    submit_dft(cfg, it)
    # Check if there are successfully converged structures, if not run the genetic algorithm again
    relaxed_file_path=Path(f"data/iter{it:03d}/dft_relaxed.xyz")


    def is_enough(path, min_structs=5):
        if not path.exists():
            return False
        structures= read(path, index=":")
        return len(structures)> min_structs

    max_retries=3
    retries=0

    while not is_enough(relaxed_file_path):
        if retries<max_retries:
            print("Number of dft relaxed structures is not enough, running the GA again!")
            run_ga(cfg,it)
            submit_dft(cfg, it)
            retries+=1
        else:
            print("Starting the GA with a different random population!")
            create_db(cfg,it)
            run_ga(cfg,it)
            submit_dft(cfg, it)
            

            # Check again — if still not enough, skip iteration completely
            if not is_enough(relaxed_file_path):
                print(f"[Iter {it}] Skipping merge_datasets due to insufficient structures.")
                break
            else:
                print(f"[Iter {it}] Success after regenerating random population.")
                break   # exit while loop
    else:
        # This runs only if the while loop exited normally (not via break)
        print(f"[Iter {it}] GA successful on retry {retries}.")


    if is_enough(relaxed_file_path):
        print("Merging the datasets...")
        merge_datasets(it, include_failed=False)
    """ print("The GA was successfull! Merging the datasets...")
    # Merge dft labels and existing data
    merge_datasets(it, include_failed=False) """


# plot mae vs iteration index
plt.figure(figsize=(6,4))
# Plot both lists
plt.plot(first_average_errors,  marker="o", label="First TEST MAE")
plt.plot(second_average_errors, marker="s", label="Second TEST MAE")
# Labels & title
plt.title("TEST MAE values across runs")
plt.xlabel("Run index")
plt.ylabel("MAE (meV/atom)")
plt.legend()

# Save to file
plt.savefig("test_mae_plot.png", dpi=300, bbox_inches="tight")
