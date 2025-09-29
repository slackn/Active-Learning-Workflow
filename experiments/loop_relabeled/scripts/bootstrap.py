#!/usr/bin/env python3
import argparse, json, random, re, yaml
from pathlib import Path
from collections import defaultdict
from typing import List, Tuple, Dict, Iterator, Optional

# -------- helpers for XYZ I/O & charge parsing --------
CHARGE_PATTERNS = [
    re.compile(r"(?:charge|q)\s*[:=]\s*([+-]?\d+)", re.I),  # charge=+1, q:-2
    re.compile(r"\bq([+-]?\d+)\b", re.I),                   # q+1, q-2
]

def iter_xyz_frames(path: Path) -> Iterator[Tuple[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        while True:
            n_line = f.readline()
            if not n_line:
                break
            n = int(n_line.strip())
            comment = f.readline().rstrip("\n")
            atoms = [f.readline().rstrip("\n") for _ in range(n)]
            frame_text = "\n".join([str(n), comment, *atoms]) + "\n"
            yield frame_text, comment

def extract_charge(comment_line: str) -> int:
    for pat in CHARGE_PATTERNS:
        m = pat.search(comment_line)
        if m:
            return int(m.group(1))
    raise ValueError(f"Could not parse charge from comment:\n{comment_line}")

def write_xyz(frames: List[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as w:
        for fr in frames:
            if not fr.endswith("\n"):
                fr += "\n"
            w.write(fr)

def stratified_split(frames_with_charge, val_frac, min_per_charge, seed):
    rng = random.Random(seed)
    buckets: Dict[int, List[str]] = defaultdict(list)
    for fr, q in frames_with_charge:
        buckets[q].append(fr)

    valid, pool, stats = [], [], {}
    for q, frs in buckets.items():
        rng.shuffle(frs)
        k = max(min_per_charge, round(len(frs) * val_frac)) if 0 < val_frac < 1 else int(val_frac)
        k = min(k, len(frs))
        valid.extend(frs[:k])
        pool.extend(frs[k:])
        stats[q] = {"total": len(frs), "valid": k, "train_pool": len(frs) - k}
    return valid, pool, stats

def make_bootstraps(train_pool, n_boot, boot_size, outdir: Path, seed):
    rng = random.Random(seed)
    N = len(train_pool)
    if N == 0:
        raise ValueError("Empty training pool after validation split.")
    size = boot_size if boot_size > 0 else N
    meta = []
    for i in range(1, n_boot + 1):
        sample = [train_pool[rng.randrange(N)] for _ in range(size)]
        out_path = outdir / f"train_boot_{i:03d}.xyz"
        write_xyz(sample, out_path)
        meta.append({"idx": i, "size": size, "path": str(out_path)})
    return meta

# -------- config-driven runner --------
def _paths_from_cfg(cfg: dict, iteration: Optional[int]) -> tuple[Path, Path]:
    ds_path = Path(cfg["data"]["dataset_pattern"].format(iter=iteration))
    out_dir = Path(cfg["data"]["iterdir_pattern"].format(iter=iteration))
    return ds_path, out_dir
    # If iteration is None, pick the latest dataset_iterXXX.xyz present.
"""     if iteration is None:
        pattern = cfg["data"]["dataset_pattern"]
        root_glob = pattern.replace("{iter:03d}", "*")
        candidates = sorted(Path().glob(root_glob))
        if not candidates:
            raise FileNotFoundError(f"No datasets found matching {root_glob}")
        # infer XXX from filenames and choose max
        def iter_from_name(p: Path) -> int:
            s = p.stem  # e.g., dataset_iter012
            digits = "".join(ch for ch in s if ch.isdigit())
            return int(digits[-3:]) if digits else -1
        latest = max(candidates, key=iter_from_name)
        # extract iteration again properly
        it = iter_from_name(latest)
        iteration = it """

    

def run_bootstrap(cfg: dict, iteration: Optional[int] = None) -> dict:
    ds_path, out_dir = _paths_from_cfg(cfg, iteration)
    prep = cfg["prep"]
    seed = int(prep.get("seed", 42)) + (iteration or 0)

    frames_with_charge = [(fr, extract_charge(c)) for fr, c in iter_xyz_frames(ds_path)]
    valid, pool, per_charge = stratified_split(
        frames_with_charge,
        val_frac=float(prep["val_frac"]),
        min_per_charge=int(prep["min_per_charge"]),
        seed=seed,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    write_xyz(valid, out_dir / "valid.xyz")
    write_xyz(pool,  out_dir / "train_pool.xyz")

    boots_meta = make_bootstraps(
        train_pool=pool,
        n_boot=int(prep["n_boot"]),
        boot_size=int(prep["boot_size"]),
        outdir=out_dir,
        seed=seed + 1,
    )

    manifest = {
        "iteration": iteration,
        "input": str(ds_path),
        "outdir": str(out_dir),
        "total_frames": len(frames_with_charge),
        "valid_frames": len(valid),
        "train_pool_frames": len(pool),
        "per_charge": per_charge,
        "n_boot": int(prep["n_boot"]),
        "boot_size": boots_meta[0]["size"] if boots_meta else 0,
        "seed_base": seed,
        "outputs": {
            "valid_xyz": str(out_dir / "valid.xyz"),
            "train_pool_xyz": str(out_dir / "train_pool.xyz"),
            "boots": boots_meta,
        },
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest

def main():
    ap = argparse.ArgumentParser(description="Config-driven stratified validation + bootstraps")
    ap.add_argument("--config", "-c", required=True, help="Path to config.yaml")
    ap.add_argument("--iter", type=int, default=None, help="Iteration index (optional). If omitted, uses latest dataset.")
    args = ap.parse_args()

    with open(args.config, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    manifest = run_bootstrap(cfg, args.iter)
    print(json.dumps(manifest, indent=2))

if __name__ == "__main__":
    main()
