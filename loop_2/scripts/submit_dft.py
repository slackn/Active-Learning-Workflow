#!/usr/bin/env python3
# submit_dft.py
import argparse
from pathlib import Path
import yaml
from ase.io import read, write
from ase.optimize import BFGS
from ase.calculators.turbomole import Turbomole
from ase.calculators.singlepoint import SinglePointCalculator
import os
from multiprocessing import Pool 



""" def get_dft_calculator():
    return Turbomole(**TM_PARAMS) """

# Turbomole parameters as a global variable
TM_PARAMS=None
charge=None
fmax=None
steps=None

def relax_one(args):
    i, a, base_iterdir= args    #get the index, atoms, and iteration folder

    folder= base_iterdir/f"atoms_{i:03d}"
    folder.mkdir(parents=True, exist_ok=True)

    # Save cwd and temporarily change to folder 
    old_cwd=Path.cwd()
    os.chdir(folder)

    confid = a.info.get("confid", "N/A")
    a.calc= Turbomole(**TM_PARAMS)
    dyn = BFGS(a, logfile=f"dft{i}_opt.log")
    try:
        print(f"[DFT] Relaxing confid={confid}")
        dyn.run(fmax=fmax, steps=steps)
        E = a.get_potential_energy()
        # keep handy metadata in XYZ comment
        a.info["energy"] = float(E)
        a.info["confid"] = confid
        #a.info["element"] = element
        a.info["charge"]  = charge
        #ok.append(a)
        ok_flag=True
        print(f"[DFT] confid={confid}, atoms_{i:03d} converged; E={E:.6f} eV")
    except RuntimeError as e:
        print(f"[DFT] confid={confid}, atoms_{i:03d} FAILED (SCF): {e}")
        # write a placeholder energy so file remains readable
        a.calc = SinglePointCalculator(a, energy=1e6)
        a.info["confid"] = confid
        #bad.append(a)
        ok_flag=False
    finally:
        os.chdir(old_cwd)
    return ok_flag, a

def submit_dft(cfg, iteration):
    global fmax, steps, TM_PARAMS, charge  #declare global variables to assign 
    n_atoms = cfg["initialization"]["n_atoms"]
    charge  = cfg["initialization"]["charge"]
    element = cfg["initialization"]["element"]

    iterdir = Path(f"data/iter{iteration:03d}")
    in_xyz  = iterdir / "selected_for_dft.xyz"
    out_ok  = iterdir / "dft_relaxed.xyz"
    out_bad = iterdir / "dft_failed.xyz"

    TM_PARAMS = {
        "total charge": cfg["dft"]["total_charge"],
        "multiplicity": cfg["dft"]["multiplicity"],
        "scf iterations": cfg["dft"]["scf_iter"],
        "basis set name": cfg["dft"]["basis_set"],
        "density functional": cfg["dft"]["density_func"],
    }

    fmax=cfg["dft"]["fmax"]
    steps=cfg["dft"]["opt_steps"]

    if not in_xyz.exists():
        raise FileNotFoundError(f"Missing input: {in_xyz}")

    frames = read(str(in_xyz), ":")
    # frame sis a python list of Atoms objects
    print(f"[DFT] Loaded {len(frames)} structures from {in_xyz}")

    # Prepare argument for pool function 
    frames_with_idx= [(i, a, iterdir) for i, a in enumerate(frames, 1)]

    ok, bad = [], []
    i=0

    # try multiprocessing
    p= Pool(10)

    results = p.map(relax_one, frames_with_idx)

    for ok_flag, result_frame in results:
        if ok_flag:
            ok.append(result_frame)
        else:
            bad.append(result_frame)



    """
    ok_flag, result_frame=p.map(relax_one, frames_with_idx)


    if ok_flag:
        ok.append(result_frame)
    elif not (ok_flag):
        bad.append(result_frame)
    """
    p.close()
    p.join()

    if ok:
        write(str(out_ok), ok)
        print(f"[DFT] Wrote {len(ok)} relaxed structures → {out_ok}")
    else:
        print("[DFT] No successful relaxations.")

    if bad:
        write(str(out_bad), bad)
        print(f"[DFT] Wrote {len(bad)} failed structures → {out_bad}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--config", required=True)
    ap.add_argument("-i", "--iter", type=int, required=True)
    args = ap.parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text())
    iteration=args.iter
    submit_dft(cfg, iteration)

if __name__ == "__main__":
    main()
