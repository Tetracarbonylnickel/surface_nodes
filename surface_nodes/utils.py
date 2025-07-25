from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import zntrack
from ipsuite import base
import ase

def mean_reduction(values, axis):
    return np.nanmean(values, axis=tuple(axis))

def max_reduction(values, axis):
    return np.nanmax(values, axis=tuple(axis))

def min_reduction(values, axis):
    return np.nanmin(values, axis=tuple(axis))

def flatten_r(values, axis):
    return values.flatten()

def check_dimension(values):
    if values.ndim > 1:
        raise ValueError(
            f"Value dimension is {values.ndim} != 1. "
            "Reduce the dimension by defining dim_reduction, "
            "use mean or max to get (n_structures,) shape."
        )

REDUCTIONS = {
    "mean": mean_reduction,
    "max": max_reduction,
    'min': min_reduction,
    'flatten': flatten_r,
}

class PlotProperties(base.IPSNode):
    data: list[ase.Atoms] = zntrack.deps()
    dim_reduction: str = zntrack.params('flatten')
    
    energy_img: Path = zntrack.outs_path(zntrack.nwd / "energy.png")
    energy_uncertainty_img: Path = zntrack.outs_path(zntrack.nwd / "energy_uncertainty.png")
    forces_img: Path = zntrack.outs_path(zntrack.nwd / "forces.png")
    forces_uncertainty_img: Path = zntrack.outs_path(zntrack.nwd / "forces_uncertainty.png")

    def get_data(self) -> list[ase.Atoms]:
        """Get the atoms data to process."""
        if self.data is not None:
            return self.data
        else:
            raise ValueError("No data given.")
        
    def run(self):
        frames = self.get_data()
        
        plot_structure_property([frame.calc.results["energy"] for frame in frames], 'energy', self.energy_img)
        plot_cartesian_property([frame.calc.results["forces"] for frame in frames], 'forces', self.dim_reduction, self.forces_img)

        if 'energy_uncertainty' in frames[0].calc.results.keys():
            plot_structure_property(
                [frame.calc.results["energy_uncertainty"] for frame in frames],
                'energy_uncertainty',
                self.energy_uncertainty_img,
            )
        else:
            plot_structure_property([frame.calc.results["energy"] for frame in frames], 'energy', self.energy_uncertainty_img)
            
        if 'energy_uncertainty' in frames[0].calc.results.keys():
            plot_cartesian_property(
                [frame.calc.results["forces_uncertainty"] for frame in frames],
                'forces_uncertainty',
                self.dim_reduction,
                self.forces_uncertainty_img,
            )
        else:
            plot_cartesian_property([frame.calc.results["forces"] for frame in frames], 'forces', self.dim_reduction, self.forces_uncertainty_img)
        
def plot_structure_property(properties, label, image):
    fig, ax = plt.subplots()
    ax.plot(properties, label=label)
    ax.set_ylabel(label)
    ax.set_xlabel("configuration")
    fig.savefig(image, bbox_inches="tight", dpi=240)


def plot_cartesian_property(properties, label, reduction, image):
    
    properties = pad_list(properties)
    
    reduction_fn = REDUCTIONS[reduction]
    properties = reduction_fn(properties, (1, 2))
    
    fig, ax = plt.subplots()
    ax.plot(properties, label=label)
    ax.set_ylabel(label)
    ax.set_xlabel("configuration")
    fig.savefig(image, bbox_inches="tight", dpi=240)
        
        
def pad_list(properties):
    
    max_rows = max(val.shape[0] for val in properties)

    # Ensure all arrays are (N, 3)
    padded_list = []
    for arr in properties:
        pad_rows = max_rows - arr.shape[0]
        if pad_rows > 0:
            # Pad with the given value along axis 0
            padding = np.full((pad_rows, 3), np.nan)
            padded_arr = np.vstack([arr, padding])
        else:
            padded_arr = arr  # Already the max size
        padded_list.append(padded_arr)
    return np.array(padded_list)
    

class SafeSelection(base.ProcessAtoms):
    
    def run(self) -> None:
        self.frames = self.get_data()

