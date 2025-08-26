#!/usr/bin/env python3
import argparse, json, subprocess, yaml
from pathlib import Path
from typing import Dict, List

def write_config(base_mace_cfg: Dict, train_file: str, valid_file: str,test_file:str, results_dir: str, name: str, out_path: Path):
    cfg = dict(base_mace_cfg)  # shallow copy is fine for flat dicts
    cfg["train_file"] = train_file
    cfg["valid_file"] = valid_file
    cfg["test_file"] = test_file 
    cfg["results_dir"] = "."
    cfg["name"] = name
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
    return out_path

def train_one(config_yaml: Path):
    # Uses mace_run_train available in your env
    subprocess.run(["mace_run_train", "--config", str(config_yaml.name)], 
    check=True,
    cwd=config_yaml.parent)

def train_ensemble_for_iteration(cfg: Dict, manifest_path: Path):
    m = json.loads(Path(manifest_path).read_text())
    it = m["iteration"]
    boots: List[Dict] = m["outputs"]["boots"]
    valid_file = m["outputs"]["valid_xyz"]

    # Fixed test set from config
    test_file = Path(cfg["data"]["test_file"]).as_posix()
    if not Path(test_file).exists():
        raise FileNotFoundError(f"Test file not found: {test_file}")

    runs_dir = Path(cfg["data"]["runsdir_pattern"].format(iter=it))
    runs_dir.mkdir(parents=True, exist_ok=True)

    base_mace_cfg = cfg["mace"]

    # Create per-boot configs
    written = []
    for b in boots:
        boot_idx = b["idx"]            # 1..N
        train_file = b["path"]         # data/iter000/train_boot_001.xyz
        model_dir = runs_dir / f"boot_{boot_idx:03d}"
        config_out = model_dir / "config.yaml"
        name = f"MACE_iter{it:03d}_boot{boot_idx:03d}"
        write_config(base_mace_cfg, train_file, valid_file,test_file, str(model_dir), name, config_out)
        written.append(config_out)

    # Sequential training (simple). You can parallelize via SLURM arrays below.
    for conf in written:
        print(f"[train] {conf}")
        train_one(conf)

def main():
    ap = argparse.ArgumentParser(description="Train MACE ensemble for a given iteration using manifest.json")
    ap.add_argument("--config", "-c", required=True, help="Path to project config.yaml")
    ap.add_argument("--manifest", required=True, help="Path to data/iterXXX/manifest.json")
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    train_ensemble_for_iteration(cfg, Path(args.manifest))

if __name__ == "__main__":
    main()
