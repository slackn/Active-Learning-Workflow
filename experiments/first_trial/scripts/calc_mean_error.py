import re
import argparse
from pathlib import Path

def extract_test_maes(logfile: Path):
    """Extract all TEST MAE E values from a MACE log file.
    Returns a list of floats (one per TEST table, in order)."""
    maes = []
    with open(logfile) as f:
        for line in f:
            if "Error-table on TEST:" in line:
                for nextline in f:
                    m = re.search(r"\|\s+\S+\s+\|\s+([\d.]+)\s+\|", nextline)
                    if m:
                        maes.append(float(m.group(1)))
                        break
    return maes  # could be [] if no TEST table


def compute_mean_test_maes_for_iteration(iter_dir: Path, iter_num: int):
    """Compute mean of first and second TEST MAE across all runs in an iteration."""
    first_maes, second_maes = [], []
    pattern = f"MACE_iter{iter_num:03d}_boot*_run-*.log"

    for logfile in iter_dir.glob(f"boot_*/logs/{pattern}"):
        if "debug" in logfile.name:
            continue  # skip debug logs
        maes = extract_test_maes(logfile)
        if len(maes) >= 1:
            first_maes.append(maes[0])
        if len(maes) >= 2:
            second_maes.append(maes[1])

    first_mean = sum(first_maes) / len(first_maes) if first_maes else None
    second_mean = sum(second_maes) / len(second_maes) if second_maes else None

    print("Test error mean: ", first_mean)
    print("Second error mean: ", second_mean)
    return first_mean, second_mean


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--iter", "-i", type=int, required=True,
                        help="Iteration number to analyze")

    args = parser.parse_args()

    iter_dir = Path("runs") / f"iter{args.iter:03d}"
    first_mean, second_mean = compute_mean_test_maes_for_iteration(iter_dir, args.iter)

    if first_mean is not None:
        print(f"[Iter {args.iter}] First TEST MAE  = {first_mean:.2f} meV/atom")
    else:
        print(f"[Iter {args.iter}] No first TEST results found.")

    if second_mean is not None:
        print(f"[Iter {args.iter}] Second TEST MAE = {second_mean:.2f} meV/atom")
    else:
        print(f"[Iter {args.iter}] No second TEST results found.")


if __name__ == "__main__":
    main()
