import typing as t

import ase
import matplotlib.pyplot as plt
import numpy as np
import tqdm
import zntrack
from ipsuite import base, models

import surface_nodes.utils as utils


class PropertyFilter(base.IPSNode):
    data: list[ase.Atoms] = zntrack.deps()
    model: models.MLModel = zntrack.deps(None)
    reference: str = zntrack.params("energy")
    cutoffs: t.Union[t.List[float]] = zntrack.params()
    direction: t.Literal["above", "below", "both"] = zntrack.params("both")
    n_configurations: int = zntrack.params(None)
    min_distance: int = zntrack.params(1)
    dim_reduction: str = zntrack.params(None)
    reduction_axis: tuple[int, ...] = zntrack.params(None)#(1, 2))

    filtered_indices: list = zntrack.outs()
    selection_plot: str = zntrack.outs_path(zntrack.nwd / "slecection.png")

    def _post_init_(self):
        if self.direction not in ["above", "below", "both"]:
            raise ValueError("'direction' should be set to 'above', 'below', or 'both'.")

        return super()._post_init_()

    def pad_list(self, inputs):
        
        if self.reduction_axis:
            max_rows = max(val.shape[0] for val in inputs)

            # Ensure all arrays are (N, 3)
            padded_list = []
            for arr in inputs:
                pad_rows = max_rows - arr.shape[0]
                if pad_rows > 0:
                    # Pad with the given value along axis 0
                    padding = np.full((pad_rows, 3), np.nan)
                    padded_arr = np.vstack([arr, padding])
                else:
                    padded_arr = arr  # Already the max size
                padded_list.append(padded_arr)
            return np.array(padded_list)
        else:
            return np.array(inputs)

    def run(self) -> t.List[int]:
        
        if self.model:
            calc = self.model.get_calculator()
            pred_atoms = calc.batch_eval(self.data)
        else:
            pred_atoms = [atoms for atoms in self.data]
            
        values = [atoms.calc.results[self.reference] for atoms in pred_atoms]
        values = self.pad_list(values)

        if self.dim_reduction is not None:
            reduction_fn = utils.REDUCTIONS[self.dim_reduction]
            values = reduction_fn(values, self.reduction_axis)

        utils.check_dimension(values)

        lower_limit, upper_limit = self.cutoffs[0], self.cutoffs[1]
        self.outlier = True

        if self.direction == "above":
            pre_selection = np.array([i for i, x in enumerate(values) if x > upper_limit])
            sorting_idx = np.argsort(values[pre_selection])[::-1]
        elif self.direction == "below":
            pre_selection = np.array([i for i, x in enumerate(values) if x < lower_limit])
            sorting_idx = np.argsort(values[pre_selection])
        else:
            pre_selection = [
                i for i, x in enumerate(values) if x < lower_limit or x > upper_limit
            ]
            if pre_selection:
                pre_selection = np.array(pre_selection)
                mean = (lower_limit + upper_limit) / 2
                dist_to_mean = abs(values[pre_selection] - mean)
                sorting_idx = np.argsort(dist_to_mean)[::-1]
            else:
                self.outlier = False
                print("no outlier")

        if self.outlier:
            self.filtered_indices = self.get_selection(pre_selection[sorting_idx])
            selection_idx = np.array(self.filtered_indices)

            values = [atoms.calc.results[self.reference] for atoms in pred_atoms]
            values = self.pad_list(values)

            if self.dim_reduction is not None:
                reduction_fn = utils.REDUCTIONS[self.dim_reduction]
                values = reduction_fn(values, self.reduction_axis)

            fig, ax = plt.subplots()
            ax.plot(values, label=self.reference)
            ax.plot(selection_idx, values[selection_idx], "x", color="red")
            ax.fill_between(
                np.arange(len(values)),
                self.cutoffs[0],
                self.cutoffs[1],
                color="black",
                alpha=0.2,
                label=f"{self.reference} +- std",
            )
            ax.set_ylabel(self.reference)
            ax.set_xlabel("configuration")

            fig.savefig(self.selection_plot, bbox_inches="tight")
        else:
            self.filtered_indices = [len(self.data) + 1]
            fig, ax = plt.subplots()
            ax.plot(1, 1, label=self.reference)
            fig.savefig(self.selection_plot, bbox_inches="tight")

    def get_selection(self, indices):
        selection = []
        for idx in indices:
            # If the value is close to any of the already selected values, skip it.
            if not selection:
                selection.append(idx)
            if not any(np.abs(idx - np.array(selection)) < self.min_distance):
                selection.append(idx)
            if len(selection) == self.n_configurations:
                break

        for id, val in enumerate(selection):
            selection[id] = int(val)
        return selection

    @property
    def frames(self) -> list[ase.Atoms]:
        return [
            self.data[i] for i in range(len(self.data)) if i not in self.filtered_indices
        ]

    @property
    def excluded_frames(self):
        if self.outlier:
            return [self.data[i] for i in self.filtered_indices]
        else:
            return []



class NoNeighborFilter(base.IPSNode):
    data: list[ase.Atoms] = zntrack.deps()
    model: models.MLModel = zntrack.deps(None)
    
    filtered_indices: list = zntrack.outs()
    
    def run(self):
        
        filtered_indices = []
        if self.model:
            idx = 0
            calc = self.model.get_calculator()

            for configuration in tqdm.tqdm(self.data, ncols=70):
                configuration: ase.Atoms
                # Run calculation
                atoms = configuration.copy()
                atoms.calc = calc
                atoms.get_potential_energy()
                if np.any(atoms.calc.results['forces_uncertainty'] < 1e-4):
                    filtered_indices.append(idx)
                idx += 1
        else:
            for idx, frame in enumerate(self.data):
                if np.any(frame.calc.results['forces_uncertainty'] < 1e-4):
                    filtered_indices.append(idx)
            
        self.filtered_indices = filtered_indices
        
        
    @property
    def frames(self) -> list[ase.Atoms]:
        return [
            self.data[i] for i in range(len(self.data)) if i not in self.filtered_indices
        ]

    @property
    def excluded_frames(self):
        if self.outlier:
            return [self.data[i] for i in self.filtered_indices]
        else:
            return []

