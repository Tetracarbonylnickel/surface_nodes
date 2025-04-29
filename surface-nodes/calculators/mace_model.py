import dataclasses
import typing

from ase.calculators.calculator import Calculator
from mace.calculators import mace_mp


@dataclasses.dataclass
class MaceCalc:
    kwargs: dict[str, typing.Any] | None = None
    device: typing.Literal["cpu", "cuda"] | None = None

    def get_calculator(self, **kwargs) -> Calculator:
        if self.kwargs is not None:
            kwargs.update(self.kwargs)
        return mace_mp(**kwargs, device=self.device)
