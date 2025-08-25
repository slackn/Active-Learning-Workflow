import argparse
from ase.io import read, write
from pathlib import Path

def merge_datasets(iter_idx: int, include_failed: bool = False, dedup_key: str = "confid"):
    iterdir = Path(f"data/iter{iter_idx:03d}")
    nextdir = Path(f"data/iter{iter_idx + 1:03d}")
    nextdir.mkdir(parents=True, exist_ok=True)

    ds_in   = iterdir / f"dataset_iter{iter_idx:03d}.xyz"
    dft_ok  = iterdir / "dft_relaxed.xyz"
    dft_bad = iterdir / "dft_failed.xyz"
    out     = nextdir / f"dataset_iter{iter_idx + 1:03d}.xyz"

    def load_frames(p: Path):
        return read(str(p), ":") if p.exists() else []

    prev_frames    = load_frames(ds_in)
    relaxed_frames = load_frames(dft_ok)
    failed_frames  = load_frames(dft_bad) if include_failed else []

    key = dedup_key
    seen, merged = set(), []
    # The order mattes because it decides which to keep
    for a in (relaxed_frames + prev_frames  + failed_frames):
        val = a.info.get(key)
        if val is None or val not in seen:
            if val is not None:
                seen.add(val)
            merged.append(a)

    if not merged:
        raise RuntimeError(f"No frames to write (checked {ds_in} and {dft_ok}).")

    write(str(out), merged)
    return out

def main():
    ap= argparse.ArgumentParser(description="Merge existing dataset with new dft labeled data.")
    ap.add_argument("--iter", "-i",type=int, required=True, help="Iteration index to find correct folders")
    ap.add_argument("--include_failed","f",action="store_true", required=False, help="Include failed dft relaxations in the merge")
    args=ap.parse_args()
    iteration=args.iter
    failed=args.include_failed
    merge_datasets(iteration, failed)

if __name__ =="__main__":
    main()
