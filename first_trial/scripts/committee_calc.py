# comittee_calc.py
import glob     #to search for file path names that match specific pattern
import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from mace.calculators import MACECalculator
from pathlib import Path


class CommitteeCalculator(Calculator):
    """
    Loads only files named like:
      checkpoints/MACE_iter{iteration:03d}_boot*_run-123.model
    """
    implemented_properties= ["energy", "forces", "sigma_E_per_atom","sigma_F_mean"]

    def __init__(self, iteration: int, use_forces=True, **kwargs):
        super().__init__(**kwargs)  #optional, it forwards optional ASE Calculator init args like label, directory
        #if you dont need these arguments, call with no kwargs
        # iteration directory relative to project root
        # look inside every boot_xxx/checkpoints/ for model files
        pattern = str("runs"/ f"iter{iteration:03d}" / "boot_*" / "checkpoints" / f"MACE_iter{iteration:03d}_boot*_run-123_stagetwo.model")
        paths = sorted(glob.glob(pattern))
        #pattern=f"checkpoints/MACE_iter{iteration:03d}_boot*_run-123.model"
        #paths= sorted(glob.glob(pattern))
        if not paths:
            raise FileNotFoundError(f"No models found with pattern {pattern}")
        self.members= [MACECalculator(model_path=p, device="cpu")for p in paths]
        self.use_forces= use_forces


    def calculate(self, atoms=None, properties=("energy","forces"),system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        Es, Fs= [], []
        for model in self.members:
            Es.append(model.get_potential_energy(atoms))
            if self.use_forces:
                Fs.append(model.get_forces(atoms))


        E_arr= np.asarray(Es)
        self.results["energy"]=float(E_arr.mean())
        sigma_E_pa= float(E_arr.std()/len(atoms))

        if self.use_forces:
            F_arr=np.stack(Fs, axis=0)
            self.results["forces"]=F_arr.mean(axis=0)
            sigma_F_mean=float(F_arr.std(axis=0).mean())
        else:
            self.results["forces"]=None
            sigma_F_mean=0

        self.results["sigma_E_pa"]=sigma_E_pa
        self.results["sigma_F_mean"]=sigma_F_mean

        """         self.results["uncertainty"] = {
        "sigma_E_per_atom": sigma_E_pa,
        "sigma_F_mean": sigma_F_mean,
        "n_models": len(self.members),
        } """

        # Save uncertainty to access in GA script
        kv= atoms.info.setdefault("key_value_pairs", {})
        kv.update({
            "sigma_E_pa": float(sigma_E_pa),
            "sigma_F_mean": float(sigma_F_mean)
        })