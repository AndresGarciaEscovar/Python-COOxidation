""" File that contains the class to define the simulation parameters."""

# ------------------------------------------------------------------------------
# Imports.
# ------------------------------------------------------------------------------

# Imports: General.
import time

from dataclasses import dataclass, field
from typing import Union

# ------------------------------------------------------------------------------
# Classes.
# ------------------------------------------------------------------------------


@dataclass
class COOxidationKMCParameters:
    """ Parameters to initiliaze the simulation.

        - length: An integer that represents the length of the lattice where the
          particles live.

        - maximum_counter: An integer, or floating, point number that represents
          the maximum time or steps that the the simulation must run for.

        - repetitions: The number of times the simulation must be run to take
          the average.

        - seed: The seed for the random number generator.
    """

    length: int = 1
    maximum_counter: Union[int, float] = 0.0
    repetitions: int = 1
    seed: int = int(time.time())
    rates: dict = field(default_factory=lambda: ({
        "Oxygen adsorption rate": 2 * 1.0,
        "Oxygen desorption rate": 2 * 1.0,
        "Oxygen diffusion rate": 2 * 1.0,
        "Carbon monoxide adsorption rate": 1 * 1.0,
        "Carbon monoxide desorption rate": 1 * 1.0,
        "Carbon monoxide diffusion rate": 2 * 1.0,
        "Langmuir-Hinshelwood reaction rate": 2 * 1.0,
        "Elay-Rideal reaction rate": 1 * 1.0
    }))
