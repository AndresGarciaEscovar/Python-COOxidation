""" File that contains the class to define the simulation parameters."""

# ------------------------------------------------------------------------------
# Imports.
# ------------------------------------------------------------------------------

# Imports: General.
import time

from dataclasses import dataclass
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

        - seed: The seed for the random number generator.
    """

    length: int = 1
    maximum_counter: Union[int, float] = 0.0
    seed: int = int(time.time())
