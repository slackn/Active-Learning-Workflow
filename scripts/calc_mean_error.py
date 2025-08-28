import re
import argparse
from pathlib import Path

def extract_first_test_mae(logfile: Path):
    """Extract the first TEST MAE E value from a MACE log file."""
    maes = []
    with open(logfile) as f:
        for line in f:
            if "Error-table on TEST:" in line:
                # Look for "MAE E / meV / atom" in the next few lines
                for nextline in f:
                    m = re.search(r"\|\s+\S+\s+\|\s+([\d.]+)\s+\|", nextline)
                    if m:
                        maes.append(float(m.group(1)))
                        break
    return maes[0] if maes else None

def compute_mean_test_mae_for_iteration(iter_dir: Path):
    """Compute mean TEST MAE across all runs in an iteration directory."""
    all_maes = []
    for logfile in iter_dir.glob("*/log.txt"):
        mae = extract_first_test_mae(logfile)
        if mae is not None:
            all_maes.append(mae)
    if not all_maes:
        return None
    return sum(all_maes) / len(all_maes)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--iter", "-i", type=int, required=True,
                        help="Iteration number to analyze")
    parser.add_argument("--runs_dir", "-r", type=str, required=True,
                        help="Base runs directory (where iter000, iter001, ... live)")
    args = parser.parse_args()

    iter_dir = Path(args.runs_dir) / f"iter{args.iter:03d}"
    mean_mae = compute_mean_test_mae_for_iteration(iter_dir)

    if mean_mae is not None:
        print(f"[Iter {args.iter}] Mean TEST MAE = {mean_mae:.2f} meV/atom")
    else:
        print(f"[Iter {args.iter}] No TEST results found.")

if __name__ == "__main__":
    main()
