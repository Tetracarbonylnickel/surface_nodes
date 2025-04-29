import dataclasses
import typing

import ase
import ase.constraints


@dataclasses.dataclass
class FixedAtomsConstraint:
    """Class to fix a layer of atoms within a MD
        simulation
    Attributes
    ----------
    indices: List[int]
        all atoms that will be fixed.
    """

    indices: typing.List[int]

    def get_constraint(self, atoms):
        return ase.constraints.FixAtoms(indices=self.indices)
