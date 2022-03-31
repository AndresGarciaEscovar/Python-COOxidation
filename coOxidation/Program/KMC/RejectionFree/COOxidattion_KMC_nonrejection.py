"""
    File that contains a procedural rejection free algorithm, that takes the
    ergodic (time-based) average.
 """

# //////////////////////////////////////////////////////////////////////////////
# Imports.
# //////////////////////////////////////////////////////////////////////////////

# General.
import numpy as np
import random


# //////////////////////////////////////////////////////////////////////////////
# Parameters.
# //////////////////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------------------
# Lattice Parameters.
# ------------------------------------------------------------------------------

particles = ["empty", "co", "o"]  # Particles of the system.
lattice = [particles[0] for _ in range(3)]   # List that represents the state of the system., i.e., contains particles.

# ------------------------------------------------------------------------------
# Random Parameters.
# ------------------------------------------------------------------------------

seed = 1234  # The seed for the random number generator.
random.seed(seed)  # Seed the random number generator.

# ------------------------------------------------------------------------------
# Rate Parameters.
# ------------------------------------------------------------------------------

# Reaction.
kLH = 1.0  # Langmuir-Hinshelwood reaction rate.
lER = 1.0  # Elay-Rideal reaction rate.

# Carbon Monoxide.
kCOdes = 1.0  # Carbon monoxide desorption rate.
kCOdif = 1.0  # Carbon monoxide diffusion rate.
kCOads = 1.0  # Carbon monoxide adsorption rate.

# Oxygen.
kOdes = 1.0  # Oxygen desorption rate.
kOdif = 1.0  # Oxygen diffusion rate.
kOads = 1.0  # Oxygen adsorption rate.

# Rates for the cells, initialize to zero.
rates = [[0.0 for _ in range(8)] for _ in range(3)]

# ------------------------------------------------------------------------------
# Statistics Parameters.
# ------------------------------------------------------------------------------

statistics = [0.0 for _ in range(3**3)]  # List that has the statistics of the system.
take_stats = False  # Flag to indicate if statistics must be taken.

# ------------------------------------------------------------------------------
# Time Parameters.
# ------------------------------------------------------------------------------

# Times.
time_equ = 100_000.0  # The time that the system must run to reach equilibrium.
time_sim = 200_000.0  # The time the simulation must run after equilibrium is reached.
time_elapsed = 0.0  # The elapsed time of the simulation.

# //////////////////////////////////////////////////////////////////////////////
# Functions.
# //////////////////////////////////////////////////////////////////////////////


# ------------------------------------------------------------------------------
# Action Methods.
# ------------------------------------------------------------------------------

def action_diffusion(sites0: tuple[int, int]) -> None:
    """
        Swaps the particles of the given two sites.

        :param sites0: The tuple that contains the sites whose particles must be
         swapped.
    """

    # --------------------------------------------------------------------------
    # Auxiliary Functions.
    # --------------------------------------------------------------------------

    def validate_swap0(sites1: tuple[int, int]) -> None:
        """
            Validates if the sites are valid for a swap.

            :param sites1: The tuple that contains the sites whose particles must be
             swapped.
        """

        # Check for errors.
        if not len(sites1) == 2:
            raise ValueError(f"The number of sites for particles to swap must be exactly two. Length = {len(sites1)}")

        elif sites1[0] == sites1[1]:
            raise ValueError(f"Trying to swap particles in the same site. Sites: {sites1}.")

        elif not all(map(lambda x: x in range(3), sites1)):
            raise ValueError(f"One or more sites are not in the lattice for a swap, {sites1}.")

        elif sites1[0] == sites1[1] + 1 or sites1[1] == sites1[0] + 1:
            raise ValueError(
                f"Sites that are swapping particles must be adjacent to each other. Values: sites0[0] = {sites1[0]}, "
                f"sites0[1] = {sites1[1]}."
            )

        # Check if the sites have valid particles.
        valid1 = (lattice[sites1[0]] == "empty" and not lattice[sites1[1]] == "empty")
        valid1 = valid1 or (lattice[sites1[1]] == "empty" and not lattice[sites1[0]] == "empty")

        if not valid1:
            raise ValueError(
                f"Trying to swap two empty, or non-empty, sites. "
                f"(sites0[0], particle0) = {(sites0[0], lattice[sites1[0]])} and "
                f"(sites0[1], particle1) = {(sites0[1], lattice[sites1[1]])}."
            )

    # --------------------------------------------------------------------------
    # Implementation.
    # --------------------------------------------------------------------------

    # Validate the swap.
    validate_swap0(sites0)

    # Make the swap.
    lattice[sites0[0]], lattice[sites0[1]] = lattice[sites0[1]], lattice[sites0[0]]


def action_adsorb(sites0: tuple[int, ...]) -> None:
    """
        Adsorbs oxygen or carbon monoxide, depending on the length of the sites
        array.

        :param sites0: The site(s) where the adsorption is going to take place.
    """

    # --------------------------------------------------------------------------
    # Auxiliary Functions.
    # --------------------------------------------------------------------------

    def validate_adsorb0(sites1: tuple[int, ...]) -> None:
        """
            Validates if the sites are valid for a swap.

            :param sites1: The tuple that contains the sites whose particles must be
             adsorbed.
        """

        # Check for errors.
        if len(sites1) not in [1, 2]:
            raise ValueError(f"The number of sites for particles to swap must be exactly one or two. Length = {len(sites1)}")


        if len(sites1) == 1:
            #TODO: Write the validation functions for adsorption of carbon monoxide.
            print()
        else:
            #TODO: Write the validation functions for adsorption of oxygen.


    # --------------------------------------------------------------------------
    # Implementation.
    # --------------------------------------------------------------------------

    # Validate the swap.
    validate_swap0(sites0)

    # Make the swap.
    lattice[sites0[0]], lattice[sites0[1]] = lattice[sites0[1]], lattice[sites0[0]]


# ------------------------------------------------------------------------------
# Get Methods.
# ------------------------------------------------------------------------------


def get_particle_id(particle0: str) -> int:
    """
        Gets the integer numerical id of the particle.

        :param particle0: The string representation of the particle.
    """

    # Look for the particle.
    for i, particle0_ in enumerate(particles):
        if particle0_ == particle0:
            return int(i)

    # Must not reach this place.
    raise ValueError(f"The id of the requested particle {particle0} cannot be found. Valid particles: {particles}.")


def get_random_number(num0_0: float, num0_1: float) -> float:
    """
        Gets a random number such that (min(num0_0,num0_1), max(num0_0,num0_1)).
        If num0 = num1, it raises an error.

        :param num0_0: The zeroth number in the range.

        :param num0_1: The first number in the range.

        :return: A random number in the range (min(num0,num1), max(num0,num1)).
    """

    # Check that the numbers are valid.
    if np.isclose(num0_0, num0_1):
        raise ValueError("The values provided for the random generator are the same. This cannot happen.")

    # Get the smallest number.
    lower0 = min(num0_0, num0_1)
    upper0 = max(num0_0, num0_1)

    # Choose a random number.
    random_number0 = random.random()
    while not 0.0 < random_number0 < 1.0:
        random_number0 = random.random()

    return lower0 + random_number0 * (upper0 - lower0)


# ------------------------------------------------------------------------------
# Record Methods.
# ------------------------------------------------------------------------------


def record_statistic(delta0: float) -> None:
    """
        Records the statistics related to the state of the lattice.

        :param delta0: The time that has gone by since the last move.
    """

    # Set the sum to zero.
    sum0 = 0

    # Get the scalar site where the statistics must be taken.
    for i, particle0 in enumerate(lattice):
        sum0 += get_particle_id(particle0) * (3 ** i)

    # Add the time to the proper state.
    statistics[sum0] += delta0


# ------------------------------------------------------------------------------
# Update Methods.
# ------------------------------------------------------------------------------


def update_rates() -> None:
    """ Updates the rates of the cells in the system. """


# //////////////////////////////////////////////////////////////////////////////
# Main Program.
# //////////////////////////////////////////////////////////////////////////////

if __name__ == "__main__":

