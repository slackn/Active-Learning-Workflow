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
    if file_path.exists():
         # Read all structures in the XYZ file
        structures = read(relaxed_file_path, index=":")  # index=":" reads all frames
        if len(structures) > 5:
            print(f"{relaxed_file_path} exists and has {len(structures)} structures → OK")
        else:
            print(f"{relaxed_file_path} exists but has only {len(structures)} structures. Running GA again!")
            run_ga(cfg, it)
    else:
        print(f"{relaxed_file_path} does not exist! Running the genetic algorithm again!")
        run_ga(cfg, it)
    # Merge dft labels and existing data
    merge_datasets(it, include_failed=False)


# plot mae vs iteration index
plt.figure(figsize=(6,4))
# Plot both lists
plt.plot(first_maes,  marker="o", label="First TEST MAE")
plt.plot(second_maes, marker="s", label="Second TEST MAE")
# Labels & title
plt.title("TEST MAE values across runs")
plt.xlabel("Run index")
plt.ylabel("MAE (meV/atom)")
plt.legend()

# Save to file
plt.savefig("test_mae_plot.png", dpi=300, bbox_inches="tight")
