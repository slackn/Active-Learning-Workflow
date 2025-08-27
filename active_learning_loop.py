import argparse
import yaml 

from scripts.bootstrap import run_bootstrap
from scripts.train_mace import train_ensemble_for_iteration
from scripts.create_db import create_db
from scripts.run_ga import run_ga
from scripts.submit_dft import submit_dft
from scripts.merge import merge_datasets
from pathlib import Path


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


for it in range(0, cfg["active_learning"]["iterations"]):
    print(it)
    if it!=0:
        #Prepare bootstrapped training datasets
        manifest=run_bootstrap(cfg, iteration=it)

        manifest_path=Path(manifest["outdir"]) / "manifest.json"
        # Train ensemble for this iteration 
        train_ensemble_for_iteration(cfg, manifest_path)
        # Create initial random population
        create_db(cfg, it)
    # Run genetic algorithm and select uncertain candidates
    run_ga(cfg, it)
    # Submit dft labeling 
    submit_dft(cfg, it)
    # Merge dft labels and existing data
    merge_datasets(it, include_failed=False)
