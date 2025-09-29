# Active Learning Workflow for MACE + GA + DFT

This repository contains an end-to-end **active learning workflow** for training MACE (Machine-learning Atomic Cluster Expansion) models with DFT-labeled datasets and candidate generation using a Genetic Algorithm (GA). It automates dataset preparation, bootstrapping, model training, uncertainty estimation, and selection of new candidates for DFT evaluation.

---

## Repository Structure
Active-Learning-Workflow/
│
├── README.md
├── requirements.txt
├── config.yaml # Main configuration file
├── scripts/ # Core scripts for workflow
│ ├── bootstrap.py # Stratified validation & bootstrap generation
│ ├── train_mace.py # Train MACE ensembles for a given iteration
│ ├── calc_mean_error.py # Compute mean TEST MAE from log files
│ ├── create_db.py # Create initial random GA population
│ ├── run_ga.py # Run GA for candidate selection
│ ├── submit_dft.py # Submit DFT relaxations
│ ├── merge.py # Merge new DFT-labeled structures into dataset
│ ├── committee_calc.py # ASE calculator for MACE ensemble
│ ├── compute_distances.py # Compute interatomic distance matrix
│ └── active_learning_workflow.py # Full active learning loop
├── data/ # Input and output datasets
└── experiments/ # Stores run-specific outputs


---

## Requirements

Install dependencies via:

```bash
pip install -r requirements.txt


Core dependencies:

numpy

PyYAML

ase

matplotlib

mace (MACE ML framework)

Configuration

All scripts use a single config.yaml file for project-wide settings. Key sections include:

data: Dataset and run directories

prep: Bootstrapping and validation parameters

initialization: GA population parameters

ga: Genetic algorithm hyperparameters

mace: MACE model training settings

active_learning: Number of active learning iterations

Usage
Full Active Learning Loop

Run the entire workflow:

python scripts/active_learning_workflow.py --config config.yaml


This will:

Generate stratified validation and bootstrap training sets.

Train MACE ensembles for each iteration.

Compute mean TEST MAE for evaluation.

Create initial random populations for GA.

Run GA to generate candidate structures.

Submit DFT relaxations.

Merge new DFT-labeled structures into the dataset.

Plot TEST MAE across iterations (test_mae_plot.png).

Individual Scripts

Bootstrap datasets:

python scripts/generate_bootstraps.py --config config.yaml --iter 0


Train MACE ensemble:

python scripts/train_mace.py --config config.yaml --manifest data/iter000/manifest.json


Compute TEST MAE:

python scripts/compute_test_mae.py --iter 0


Create initial GA database:

python scripts/create_db.py --config config.yaml --iter 0


Run GA:

python scripts/run_ga.py --config config.yaml --iter 0


Submit DFT relaxations:

python scripts/submit_dft.py --config config.yaml --iter 0


Merge datasets:

python scripts/merge.py --iter 0

Workflow Overview

Bootstrapping: Create stratified validation and training subsets per charge.

Training: Train ensemble MACE models for each bootstrap.

Uncertainty Estimation: Use CommitteeCalculator to compute ensemble-averaged energies and forces with uncertainties.

GA Candidate Generation: Propose new structures using genetic operations.

DFT Labeling: Relax candidate structures using DFT.

Dataset Update: Merge successfully relaxed structures into the dataset.

Iteration: Repeat for the number of active learning iterations defined in config.yaml.

Outputs

data/iterXXX/valid.xyz – Validation set

data/iterXXX/train_pool.xyz – Training pool

data/iterXXX/train_boot_XXX.xyz – Bootstrapped training sets

data/iterXXX/*.db – ASE databases for GA

data/iterXXX/dft_relaxed.xyz – Successfully relaxed DFT structures

manifest.json – Metadata for each iteration

test_mae_plot.png – TEST MAE vs iteration index

Notes


Ensure mace_run_train is available in your environment PATH.
