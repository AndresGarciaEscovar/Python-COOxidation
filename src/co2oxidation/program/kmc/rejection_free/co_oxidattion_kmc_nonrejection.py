"""
    File that contains a procedural rejection free algorithm, that takes the
    ergodic (time-based) average.
 """

# //////////////////////////////////////////////////////////////////////////////
# Imports.
# //////////////////////////////////////////////////////////////////////////////


# General.
import copy as cp
import csv
import os
import random

from datetime import datetime
from itertools import product

# Third party.
import numpy as np

from matplotlib import pyplot as plt


# //////////////////////////////////////////////////////////////////////////////
# Parameters.
# //////////////////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------------------
# File Parameters.
# ------------------------------------------------------------------------------

folder_name = os.path.dirname(__file__) + os.sep
file_name = "COOxidation"

# ------------------------------------------------------------------------------
# Lattice Parameters.
# ------------------------------------------------------------------------------

# Physical parameters.
area = 2.5 * 10**-19  # Area of the substrate (m^2).
temperature = 600  # Temperature of the heat bath (K).
factor_arrhenius = 10**13  # Arrhenius rate pre-factor (s^-1).

# For the simulation.
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
eact_lh = 1.0  # Langmuir-Hinshelwood activation energy (ev).

# Carbon Monoxide.
eact_co_des = 1.0  # Carbon monoxide activation energy (ev).
eact_co_dif = 1.0  # Carbon monoxide diffusion energy (ev).
eact_co_er = 1.0  # Elay-Rideal activation energy for carbon monoxide (ev).

partial_pressure_co = 10.0**5  # Partial Pressure for carbon monoxide (Pa).
molar_mass_co = 0.016 + 0.012  # Molar mass of carbon monoxide (CO) (kg/mol).

# Oxygen.
eact_o_des = 1.0  # Oxygen desorption activation energy (ev).
eact_o_dif = 1.0  # Oxygen diffusion activation energy (ev).
eact_o_er = 1.0  # Elay-Rideal activation energy for oxygen (ev).

partial_pressure_o = 10.0**5  # Partial Pressure for oxygen (Pa).
molar_mass_o = 0.016 + 0.016  # Molar mass of oxygen molecule (O2) (kg/mol).

# Rates for the cells, initialize to zero.
rates = [[0.0 for _ in range(9)] for _ in range(3)]

# ------------------------------------------------------------------------------
# Statistics Parameters.
# ------------------------------------------------------------------------------

statistics = [0.0 for _ in range(3**3)]  # List that has the statistics of the system.
take_stats = False  # Flag to indicate if statistics must be taken.

# ------------------------------------------------------------------------------
# Time Parameters.
# ------------------------------------------------------------------------------

# Times.
time_equ = 1_000_000.0  # The time that the system must run to reach equilibrium.
time_sim = 2_000_000.0  # The time the simulation must run after equilibrium is reached.
time_elapsed = 0.0  # The elapsed time of the simulation.

# //////////////////////////////////////////////////////////////////////////////
# Functions.
# //////////////////////////////////////////////////////////////////////////////


# ------------------------------------------------------------------------------
# Action Methods.
# ------------------------------------------------------------------------------

def action_adsorption(site0: int, particle0: str) -> None:
    """
        Adsorbs oxygen or carbon monoxide, using the given site as reference.

        :param site0: The site where the adsorption is going to take place.

        :param particle0: The type of particle to be adsorbed; either oxygen 'o'
         or carbon monoxide 'co'.
    """

    # --------------------------------------------------------------------------
    # Auxiliary Functions.
    # --------------------------------------------------------------------------

    def validate_adsorption0(site1: int, particle1: str) -> None:
        """
            Validates if the site is valid for an adsorption of the given
            particle.

            :param site1: The site where the adsorption is going to take place.

            :param particle1: The type of particle to be adsorbed; either oxygen
             'o' or carbon monoxide 'co'.
        """

        # Check for type errors.
        if type(site1) is not int:
            raise ValueError(f"The site must be an integer. Current type: {type(site1)}.")

        # Validate the particle type.
        if particle1 not in ["o", "co"]:
            raise ValueError(
                f"Particle is not valid for adsorption. Must be an oxygen (o) or carbon monoxide (co). Current "
                f"particle: {particle1}"
            )

        # Auxiliary variables.
        range1 = list(range(2)) if particle1 == "o" else list(range(3))
        site1_ = site1 + 1 if particle1 == "o" else site1
        sites_adsorption1 = [site1, site1_] if particle1 == "o" else [site1]

        # Check the adsorption sites and lattice configuration at the sites.
        if site1 not in range1:
            raise ValueError(
                f"The site for adsorption of '{particle1}' must be in the range of integers {range1}. Requested site: "
                f"{site1}."
            )

        elif not all(map(lambda x: lattice[x] == "empty", sites_adsorption1)):
            raise ValueError(
                f"The site for adsorption of '{particle1}' must be in the range of integers {range1} and the sites "
                f"must be empty. Requested site: {site1}, Current site particle(s) {lattice[site1: site1_ + 1]}."
            )

    # --------------------------------------------------------------------------
    # Implementation.
    # --------------------------------------------------------------------------

    # Validate the swap.
    validate_adsorption0(site0, particle0)

    # Adsorb the corresponding particles.
    if particle0 == "o":
        lattice[site0], lattice[site0 + 1] = "o", "o"
        return
    lattice[site0] = "co"


def action_desorption(site0: int, particle0: str) -> None:
    """
        Desorbs oxygen or carbon monoxide, using the given site as reference.

        :param site0: The site where the desorption is going to take place.

        :param particle0: The type of particle to be desorbed; either oxygen 'o'
         or carbon monoxide 'co'.
    """

    # --------------------------------------------------------------------------
    # Auxiliary Functions.
    # --------------------------------------------------------------------------

    def validate_desorption0(site1: int, particle1: str) -> None:
        """
            Validates if the site is valid for a desorption of the given
            particle.

            :param site1: The site where the desorption is going to take place.

            :param particle1: The type of particle to be desorbed; either oxygen
             'o' or carbon monoxide 'co'.
        """

        # Check for type errors.
        if type(site1) is not int:
            raise ValueError(f"The site must be an integer. Current type: {type(site1)}.")

        # Validate the particle type.
        if particle1 not in ["o", "co"]:
            raise ValueError(
                f"Particle is not valid for adsorption. Must be an oxygen (o) or carbon monoxide (co). Current "
                f"particle: {particle1}"
            )

        # Auxiliary variables.
        range1 = list(range(2)) if particle1 == "o" else list(range(3))
        site1_ = site1 + 1 if particle1 == "o" else site1
        sites_adsorption1 = [site1, site1_] if particle1 == "o" else [site1]

        # Check the adsorption sites and lattice configuration at the sites.
        if site1 not in range1:
            raise ValueError(
                f"The site for adsorption of '{particle1}' must be in the range of integers {range1}. Requested site: "
                f"{site1}."
            )

        elif particle1 == 'o' and not all(map(lambda x: lattice[x] == "o", sites_adsorption1)):
            raise ValueError(
                f"The site for desorption of '{particle1}' must be in the range of integers {range1} and the sites "
                f"must be oxygen. Requested site: {site1}, Current site particle(s) {lattice[site1: site1_ + 1]}."
            )

        elif particle1 == 'co' and not lattice[site1] == "co":
            raise ValueError(
                f"The site for desorption of '{particle1}' must be in the range of integers {range1} and the sites "
                f"must be oxygen. Requested site: {site1}, Current site particle(s) {lattice[site1]}."
            )

    # --------------------------------------------------------------------------
    # Implementation.
    # --------------------------------------------------------------------------

    # Validate the swap.
    validate_desorption0(site0, particle0)

    # Desorb the corresponding particles.
    if particle0 == "o":
        lattice[site0], lattice[site0 + 1] = "empty", "empty"
        return
    lattice[site0] = "empty"


def action_diffusion(site0: int) -> None:
    """
        Swaps the particles of the given site with the one adjacent to the right
        of it.

        :param site0: The numerical identification of the site where the
         particles must swap.
    """

    # --------------------------------------------------------------------------
    # Auxiliary Functions.
    # --------------------------------------------------------------------------

    def validate_diffusion0(site1: int) -> None:
        """
            Validates if the site is valid for a swap.

            :param site1: The numerical identification of the site where the
             particles must swap.
        """

        # Check for errors.
        if site1 not in range(2) or type(site1) is not int:
            raise ValueError(f"The site must be an integer in the range {list(range(2))}; requested site: {site1}.")

        # Check if the sites have valid particles to swap.
        valid1 = (lattice[site1] == "empty" and not lattice[site1 + 1] == "empty")
        valid1 = valid1 or (lattice[site1 + 1] == "empty" and not lattice[site1] == "empty")

        if not valid1:
            raise ValueError(
                f"Trying to swap two empty, or non-empty, sites. (site, particle) = {(site1, lattice[site1])} and "
                f"(site + 1, particle) = {(site1 + 1, lattice[site1 + 1])}."
            )

    # --------------------------------------------------------------------------
    # Implementation.
    # --------------------------------------------------------------------------

    # Validate the swap.
    validate_diffusion0(site0)

    # Make the swap.
    lattice[site0], lattice[site0 + 1] = lattice[site0 + 1], lattice[site0]


def action_reaction_er(site0: int, particle0: str = "o") -> None:
    """
        Performs a standard Elay-Rideal reaction of oxygen in the system; i.e.,
        an oxygen particle on the substrate associatively desorbs with a carbon
        monoxide in the environment.

        :param site0: The site where the reaction is going to take place.

        :param particle0: The particle to be desorbed through the Elay-Rideal
         process.
    """

    # --------------------------------------------------------------------------
    # Auxiliary Functions.
    # --------------------------------------------------------------------------

    def validate_reaction_er0(site1: int, particle1: str) -> None:
        """
            Validates if the site is valid for a standard Elay-Rideal reaction
            of oxygen in the system; i.e., an oxygen particle on the substrate
            associatively desorbs with a carbon monoxide in the environment.

            :param site1: The site where the reaction is going to take place.

            :param particle1: The particle to be desorbed through the
             Elay-Rideal process.
        """

        # Check for type errors.
        if type(site1) is not int:
            raise ValueError(f"The site must be an integer. Current type: {type(site1)}.")

        # Check the adsorption sites and lattice configuration at the sites.
        if site1 not in list(range(3)):
            raise ValueError(
                f"The site for the Elay-Rideal reaction must be in the range of integers {list(range(2))}. "
                f"Requested site: {site1}."
            )

        # Check the site configuration.
        if not lattice[site1] == particle1:
            raise ValueError(
                f"The site in the Elay-Rideal reaction must have an oxygen. Requested site: {site1}, Current site "
                f"particle(s) {lattice[site1]}."
            )

    # --------------------------------------------------------------------------
    # Implementation.
    # --------------------------------------------------------------------------

    # Validate the swap.
    validate_reaction_er0(site0, particle0)

    # Perform the reaction.
    lattice[site0] = "empty"


def action_reaction_lh(site0: int) -> None:
    """
        Performs a standard Langmuir-Hinshelwood reaction in the system; i.e.,
        two neighboring particles (o-cO or co-o) on the substrate associatively
        desorb.

        :param site0: The left-most site where the reaction is going to take
         place.
    """

    # --------------------------------------------------------------------------
    # Auxiliary Functions.
    # --------------------------------------------------------------------------

    def validate_reaction_lh0(site1: int) -> None:
        """
            Validates if the site is valid for a standard Langmuir-Hinshelwood
            reaction in the system; i.e., two neighboring particles
            (o-cO or co-o) associatively desorb.

            :param site1: The left-most site where the reaction is going to take
             place.
        """

        # Check for type errors.
        if type(site1) is not int:
            raise ValueError(f"The site must be an integer. Current type: {type(site1)}.")

        # Check the adsorption sites and lattice configuration at the sites.
        if site1 not in list(range(2)):
            raise ValueError(
                f"The site for the Langmuir-Hinshelwood reaction must be in the range of integers {list(range(2))}. "
                f"Requested site: {site1}."
            )

        # Check the site configuration.
        valid1 = lattice[site1] == "o" and lattice[site1 + 1] == "co"
        valid1 = valid1 or (lattice[site1 + 1] == "o" and lattice[site1] == "co")

        if not valid1:
            raise ValueError(
                f"The sites in the Langmuir-Hinshelwood reaction must have an oxygen and a carbon monoxide, in any "
                f"order. Requested site: {site1}, Current site particle(s) {lattice[site1: site1 + 2]}."
            )

    # --------------------------------------------------------------------------
    # Implementation.
    # --------------------------------------------------------------------------

    # Validate the swap.
    validate_reaction_lh0(site0)

    # Perform the reaction.
    lattice[site0], lattice[site0 + 1] = "empty", "empty"


# ------------------------------------------------------------------------------
# Get Methods.
# ------------------------------------------------------------------------------


def get_boltzmann_constant(units0) -> float:
    """
        Gets the Boltzmann constant in the requested units. Currently, only in
        electron volts (ev) and joules (j); not case-sensitive.

        :param units0: The units in which the Boltzmann constant must be
         provided.

        :return: The value of the Boltzmann constant in the requested units.
    """

    # --------------------------------------------------------------------------
    # Auxiliary functions.
    # --------------------------------------------------------------------------

    def validate_units0(units1: str) -> None:
        """
            Validates that the units in which the constant is requested are
            valid.

            :param units1: The units in which the Boltzmann constant must be
             provided.
        """

        # Possible units.
        valid1 = ["ev", "j"]

        if units1.strip().lower() not in valid1:
            raise ValueError(f"The requested units {units1} are not valid.")

    # --------------------------------------------------------------------------
    # Implementation.
    # --------------------------------------------------------------------------

    # Validate the units.
    validate_units0(units0)

    # Return the value of the constant.
    kb0 = 8.617_333 * 10**-5 if units0.strip().lower() == "ev" else 1.380_649 * 10**-23

    return float(kb0)


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


def get_rate_adsorption(pressure0: float, area0: float, molar_mass0: float, temperature0: float = 300) -> float:
    """
        Gets the adsorption rate of the specie, given the system variables and
        parameters. The values must be provided in SI units.

        :param pressure0: The partial pressure of the gas in the environment, in
         SI units, i.e., Pascals (Pa).

        :param area0: The effective adsorption are of the substrate, in SI
         units, i.e., meters squared (m^2).

        :param molar_mass0: The molar mass of the particules in the gas, in SI
         units, i.e., kilograms per mole (kg/mol).

        :param temperature0: The temperature of the heat bath, in Kelvin.

        :return: The rate of adsorption.
    """

    # Get the Boltzmann constant in joules.
    units0 = "j"
    kb0 = get_boltzmann_constant(units0)

    # Calculate the denominator.
    power0 = 0.5
    denominator0 = (2.0 * np.pi * molar_mass0 * kb0 * temperature0)**power0

    return (area0 * pressure0) / denominator0


def get_rate_arrhenius(energy0: float, temperature0: float = 300, factor0: float = 10**12, units0: str = "ev") -> float:
    """
        Gets the Arrhenius rate for the process, given the activation energy.
        The parameters of the equation can also be set, as wells as the units
        for the Boltzmann constant.

        :param energy0: The activation energy of the process.

        :param temperature0: The temperature of the heat bath, in Kelvin.

        :param factor0: The pre-factor of the Arrhenius rate exponential.

        :param units0: Only options, as of now, are electron volts (ev) and
         joules (j); not case-sensitive.

        :return: The Arrhenius rate of the process, given the parameters.
    """

    # Get the exponential coefficient and Boltzmann constant.
    kb0 = get_boltzmann_constant(units0)
    coefficient0 = energy0 / (kb0 * temperature0)

    return float(factor0 * np.exp(-coefficient0))


def get_statistics(configuration0: list) -> float:
    """
        Gets the statistic of the requested configuration.

        :param configuration0: The list of particles with the requested
         configuration.

        :return: The statistic related to the given configuration.
    """

    # --------------------------------------------------------------------------
    # Auxiliary Functions.
    # --------------------------------------------------------------------------

    def validate_configuration0(configuration1: list) -> None:
        """
            Validate that the requested configuration exists.

            :param configuration1: The list of particles with the requested
             configuration.
        """

        # Check that the length is consistent.
        if not len(configuration1) == len(lattice):
            raise ValueError(
                f"The requested configuration must have the same length as that of the lattice. Lattice length: "
                f"{len(lattice)}, requested configuration length: {len(configuration1)}."
            )

        elif not all(map(lambda x: x.strip().lower() in particles, configuration1)):
            raise ValueError(
                f"The particles in the configuration are NOT valid. Configuration {configuration1}, valid particles: "
                f"{particles}."
            )

    # --------------------------------------------------------------------------
    # Implementation.
    # --------------------------------------------------------------------------

    # Validate the configuration
    validate_configuration0(configuration0)

    # Get the scalar identifier.
    sum0 = 0
    for i0, particle0 in enumerate(configuration0):
        sum0 += get_particle_id(particle0.strip().lower()) * 3**i0

    return statistics[sum0]


def get_total_rate() -> float:
    """
        Gets the total rate of the system.

        :return: The total rate of the system, i.e., the sum of the cumulative
         rates of each cell in the system.
    """

    # Total rates of the system.
    total_rate0 = sum([cell_rates0[-1] for cell_rates0 in rates])

    return total_rate0


# ------------------------------------------------------------------------------
# Plot Methods.
# ------------------------------------------------------------------------------


def plot_pies() -> None:
    """
        Plots the pie charts of the results.
    """

    # Get the information.
    information0 = [
        [
            sum([get_statistics(["empty", p0_0, p0_1]) for p0_0, p0_1 in product(particles, particles)]),
            sum([get_statistics(["co", p0_0, p0_1]) for p0_0, p0_1 in product(particles, particles)]),
            sum([get_statistics(["o", p0_0, p0_1]) for p0_0, p0_1 in product(particles, particles)]),
        ],
        [
            sum([get_statistics([p0_0, "empty", p0_1]) for p0_0, p0_1 in product(particles, particles)]),
            sum([get_statistics([p0_0, "co", p0_1]) for p0_0, p0_1 in product(particles, particles)]),
            sum([get_statistics([p0_0, "o", p0_1]) for p0_0, p0_1 in product(particles, particles)]),
        ],
        [
            sum([get_statistics([p0_0, p0_1, "empty"]) for p0_0, p0_1 in product(particles, particles)]),
            sum([get_statistics([p0_0, p0_1, "co"]) for p0_0, p0_1 in product(particles, particles)]),
            sum([get_statistics([p0_0, p0_1, "o"]) for p0_0, p0_1 in product(particles, particles)]),
        ]
    ]

    # Configure the plots.
    figure0, axes0 = plt.subplots(nrows=1, ncols=3)

    # Configure the plots.
    axes0[0].pie(information0[0], autopct='%1.1f%%')
    axes0[0].axis("equal")

    axes0[1].pie(information0[1], autopct='%1.1f%%')
    axes0[1].axis("equal")

    patches, *texts = axes0[2].pie(information0[2], autopct='%1.1f%%')
    axes0[2].axis("equal")

    figure0.legend(patches, ["Empty", "CO", "O"])
    plt.show()
    plt.clf()
    plt.close()

# ------------------------------------------------------------------------------
# Rate Methods.
# ------------------------------------------------------------------------------


def rate_adsorption(site0: int, particle0: str) -> float:
    """
        Returns the adsorption rate of the given site.

        :param site0: The site whose rate is to be calculated.

        :param particle0: The particle whose rate is to be calculated.

        :return: The adsorption rate for the given site.
    """

    # --------------------------------------------------------------------------
    # Auxiliary Functions.
    # --------------------------------------------------------------------------

    def validate_adsorption0(particle1: str) -> None:
        """
            Validates if the provided particle is valid for adsorption.

            :param particle1: The particle that is adsorbing.
        """

        # Raise error if particle doesn't exist.
        if particle1 == "empty" or particle1 not in particles:
            raise ValueError(
                f"Trying to adsorb a non-valid particle. Current particle: {particle1}, valid particles: "
                f"{particles[1:]}"
            )

    # --------------------------------------------------------------------------
    # Implementation.
    # --------------------------------------------------------------------------

    # Validate the adsorption process.
    validate_adsorption0(particle0)

    # Out of range sites have a zero rate.
    if site0 > 1 and particle0 == "o":
        return 0.0

    # Check lattice configuration.
    if particle0 == "o" and lattice[site0] == "empty" and lattice[site0 + 1] == "empty":
        partial_pressure0 = partial_pressure_o
        molar_mass0 = molar_mass_o

    elif particle0 == "co" and lattice[site0] == "empty":
        partial_pressure0 = partial_pressure_co
        molar_mass0 = molar_mass_co

    else:
        return 0.0

    return get_rate_adsorption(
        pressure0=partial_pressure0, area0=area, molar_mass0=molar_mass0, temperature0=temperature
    )


def rate_desorption(site0: int, particle0: str) -> float:
    """
        Returns the desorption rate of the particles at the given site.

        :param site0: The site whose rate is to be calculated.

        :param particle0: The particle whose rate is to be calculated.

        :return: The desorption rate for the given site.
    """

    # --------------------------------------------------------------------------
    # Auxiliary Functions.
    # --------------------------------------------------------------------------

    def validate_desorption0(particle1: str) -> None:
        """
            Validates if the provided particle is valid for desorption.

            :param particle1: The particle that is desorbing.
        """

        # Raise error if particle doesn't exist.
        if particle1 == "empty" or particle1 not in particles:
            raise ValueError(
                f"Trying to swap a non-valid particle. Current particle: {particle1}, valid particles: {particles[1:]}"
            )

    # --------------------------------------------------------------------------
    # Implementation.
    # --------------------------------------------------------------------------

    # Validate the desorption process.
    validate_desorption0(particle0)

    # Out of range sites have a zero rate.
    if site0 > 1 and particle0 == "o":
        return 0.0

    # Check lattice configuration.
    if particle0 == "o" and lattice[site0] == "o" and lattice[site0 + 1] == "o":
        eact0 = eact_o_des

    elif particle0 == "co" and lattice[site0] == "co":
        eact0 = eact_co_des

    else:
        return 0.0

    return get_rate_arrhenius(energy0=eact0, temperature0=temperature, factor0=factor_arrhenius, units0="ev")


def rate_diffusion(site0: int, particle0: str) -> float:
    """
        Returns the diffusion rate of the given particle at the given site.

        :param site0: The site whose rate is to be calculated.

        :param particle0: The particle whose rate is to be calculated.

        :return: The diffusion rate for the given particle at the given site.
    """

    # --------------------------------------------------------------------------
    # Auxiliary Functions.
    # --------------------------------------------------------------------------

    def validate_diffusion0(particle1: str) -> None:
        """
            Validates if the provided particle is valid for diffusion.

            :param particle1: The particle that is diffusing.
        """

        # Raise error if particle doesn't exist.
        if particle1 == "empty" or particle1 not in particles:
            raise ValueError(
                f"Trying to swap a non-valid particle. Current particle: {particle1}, valid particles: {particles[1:]}."
            )

    # --------------------------------------------------------------------------
    # Implementation.
    # --------------------------------------------------------------------------

    # Validate the particle.
    validate_diffusion0(particle0)

    # Out of range sites have a zero rate.
    if site0 > 1:
        return 0.0

    # Check lattice configuration.
    valid1 = lattice[site0] == "empty" and lattice[site0 + 1] == particle0
    valid1 = valid1 or (lattice[site0] == particle0 and lattice[site0 + 1] == "empty")

    if not valid1:
        return 0.0

    # Get the proper activation energy.
    valid1 = lattice[site0] == "o" or lattice[site0 + 1] == "o"
    eact0 = eact_o_dif if valid1 else eact_co_dif

    return get_rate_arrhenius(energy0=eact0, temperature0=temperature, factor0=factor_arrhenius, units0="ev")


def rate_er(site0: int, particle0: str) -> float:
    """
        Returns the Elay-Rideal reaction rate of the given site.

        :param site0: The site whose rate is to be calculated.

        :param particle0: The particle whose Elay-Rideal reaction rate is to be
         calculated.

        :return: The Elay-Rideal rate for the given site.
    """

    # --------------------------------------------------------------------------
    # Auxiliary Functions.
    # --------------------------------------------------------------------------

    def validate_er0(particle1: str) -> None:
        """
            Validates if the provided particle is valid for an Elay-Rideal
            reaction.

            :param particle1: The particle that is adsorbing.
        """

        # Raise error if particle doesn't exist.
        if particle1 == "empty" or particle1 not in particles:
            raise ValueError(
                f"Trying to perform an Elay-Rideal reaction fro p a non-valid particle. Current particle: "
                f"{particle1[1:]}, valid particles: {particles[1:]}"
            )

    # --------------------------------------------------------------------------
    # Implementation.
    # --------------------------------------------------------------------------

    # Validate to calculate the er reaction.
    validate_er0(particle0)

    # Out of range sites have a zero rate.
    if particle0 == "co" or not lattice[site0] == particle0:
        return 0.0

    # Choose the proper activation energy.
    eact0 = eact_o_er if particle0 == "o" else eact_co_er

    return get_rate_arrhenius(energy0=eact0, temperature0=temperature, factor0=factor_arrhenius, units0="ev")


def rate_lh(site0: int) -> float:
    """
        Returns the Langmuir-Hinshelwood reaction rate of the given site.

        :param site0: The site whose rate is to be calculated.

        :return: The Langmuir-Hinshelwood rate for the given site.
    """

    # Out of range sites have a zero rate.
    if site0 > 1:
        return 0.0

    # Check lattice configuration.
    valid1 = lattice[site0] == "o" and lattice[site0 + 1] == "co"
    valid1 = valid1 or (lattice[site0] == "co" and lattice[site0 + 1] == "o")

    if not valid1:
        return 0.0

    return get_rate_arrhenius(energy0=eact_lh, temperature0=temperature, factor0=factor_arrhenius, units0="ev")


# ------------------------------------------------------------------------------
# Record Methods.
# ------------------------------------------------------------------------------


def record_statistics(delta0: float) -> None:
    """
        Records the statistics related to the state of the lattice.

        :param delta0: The time that has gone by since the last move.
    """

    # Get the scalar site where the statistics must be taken.
    sum0 = 0
    for i, particle0 in enumerate(lattice):
        sum0 += get_particle_id(particle0) * (3 ** i)

    # Add the time to the proper state.
    statistics[sum0] += delta0


def record_statistics_final() -> None:
    """
        Normalizes the statistics of the simulation.
    """

    # Get the average statistics.
    normalizing0 = sum(statistics)
    for i0 in range(len(statistics)):
        statistics[i0] /= normalizing0


# ------------------------------------------------------------------------------
# Run Methods.
# ------------------------------------------------------------------------------


def run_simulation(equilibrate0: bool = False) -> tuple:
    """
        Runs the simulation of the system. If the simulation is run with the
        equilibration flag set to true, it runs the simulation without taking
        any statistics.

        :param equilibrate0: Determines if statistics must be taken while the
         simulation is running.

        :return: A tuple with the total physical elapsed time and the
         corresponding number of Monte Carlo simulation steps.
    """

    # Set the elapsed time to zero.
    t0 = 0.0
    steps0 = 0
    max_time0 = time_equ if equilibrate0 else time_sim

    # Run until the time reaches equilibrium.
    while t0 < max_time0:
        # Advance the number of steps.
        steps0 += 1

        # Stochastic time advance.
        update_rates()
        delta_t0 = -np.log(get_random_number(0.0, 1.0)) / get_total_rate()
        t0 += delta_t0

        # Record the statistics if needed.
        if not equilibrate0:
            record_statistics(delta_t0)

        # Make a move.
        cell0 = select_cell()
        process0 = select_process(cell0)
        update_lattice(cell0, process0)

    # Record the final statistics.
    if not equilibrate0:
        t0 = sum(statistics)
        record_statistics_final()

    return t0, steps0


# ------------------------------------------------------------------------------
# Save Methods.
# ------------------------------------------------------------------------------


def save_statistics(elapsed_time0: float, steps0: int, file_name0: str) -> None:
    """
        Saves the statistics to the given file.

        :param elapsed_time0: The elapsed time of the simulation after
         equilibration.

        :param steps0: The number of kinetic Monte Carlo steps that were used
         for the simulation to run.

        :param file_name0: The name of the file where the simulation is to be
         saved; without any extension.
    """

    # Date and time.
    date0 = datetime.now().strftime("%Y%m%d-%H:%M:%S")

    # Carbon monoxide parameters.
    co_adsorption_rate0 = get_rate_adsorption(
        pressure0=partial_pressure_co, area0=area, molar_mass0=molar_mass_co, temperature0=temperature
    )

    co_desorption_rate0 = get_rate_arrhenius(
        energy0=eact_co_des, temperature0=temperature, factor0=factor_arrhenius, units0="ev"
    )

    co_diffusion_rate0 = get_rate_arrhenius(
        energy0=eact_co_dif, temperature0=temperature, factor0=factor_arrhenius, units0="ev"
    )

    # Oxygen parameters.
    o_adsorption_rate0 = get_rate_adsorption(
        pressure0=partial_pressure_o, area0=area, molar_mass0=molar_mass_o, temperature0=temperature
    )

    o_desorption_rate0 = get_rate_arrhenius(
        energy0=eact_o_des, temperature0=temperature, factor0=factor_arrhenius, units0="ev"
    )

    o_diffusion_rate0 = get_rate_arrhenius(
        energy0=eact_o_dif, temperature0=temperature, factor0=factor_arrhenius, units0="ev"
    )

    # Reactions.
    lh_rate0 = get_rate_arrhenius(
        energy0=eact_lh, temperature0=temperature, factor0=factor_arrhenius, units0="ev"
    )

    er_o_reaction_rate0 = get_rate_arrhenius(
        energy0=eact_o_er, temperature0=temperature, factor0=factor_arrhenius, units0="ev"
    )

    er_co_reaction_rate0 = 0.0

    # Simulation information.
    information0 = [
        [
            "CO Oxidation Non-rejection Ergodic Average",
            f"Date: {date0}",
            f"Random generator seed: {seed}",
            f"Equilibration time (s): {time_equ}",
            f"Total elapsed time (s): {elapsed_time0}",
            f"Equivalent Monte Carlo steps: {steps0}",
            f"Simulation time (s): {time_sim}",
            f"Temperature (K): {temperature}",
            f"Surface area (m^2): {area}",
            f"Arrhenius factor (s^-1): {factor_arrhenius:.5e}",
            f"Molar mass carbon monoxide (kg/mol): {molar_mass_co}",
            f"Partial pressure carbon monoxide (Pa): {partial_pressure_co}",
            f"Molar mass oxygen (kg/mol): {molar_mass_o}",
            f"Partial pressure oxygen (Pa): {partial_pressure_o}",
            f"Langmuir-Hinshelwood reaction rate: {lh_rate0:.5e}",
            f"Elay-Rideal carbon monoxide reaction rate: {er_co_reaction_rate0:.5e}",
            f"CO adsorption rate (s^-1): {co_adsorption_rate0:.5e}",
            f"CO desorption rate (s^-1): {co_desorption_rate0:.5e}",
            f"CO diffusion rate (s^-1): {co_diffusion_rate0:.5e}",
            f"Elay-Rideal oxygen reaction rate: {er_o_reaction_rate0:.5e}",
            f"O adsorption rate (s^-1): {o_adsorption_rate0:.5e}",
            f"O desorption rate (s^-1): {o_desorption_rate0:.5e}",
            f"O diffusion rate (s^-1): {o_diffusion_rate0:.5e}"
        ]
    ]

    # Get the single site statistics.
    table_header0 = [
        "Site Number", *[
            f"{particle0.title()}" if i0 == 0 else f"{particle0.upper()}" for i0, particle0 in enumerate(particles)
        ]
    ]
    information0.append(table_header0)

    information0.extend(
        [
            [
                1,
                sum([get_statistics(["empty", p0_0, p0_1]) for p0_0, p0_1 in product(particles, particles)]),
                sum([get_statistics(["co", p0_0, p0_1]) for p0_0, p0_1 in product(particles, particles)]),
                sum([get_statistics(["o", p0_0, p0_1]) for p0_0, p0_1 in product(particles, particles)]),
            ],
            [
                2,
                sum([get_statistics([p0_0, "empty", p0_1]) for p0_0, p0_1 in product(particles, particles)]),
                sum([get_statistics([p0_0, "co", p0_1]) for p0_0, p0_1 in product(particles, particles)]),
                sum([get_statistics([p0_0, "o", p0_1]) for p0_0, p0_1 in product(particles, particles)]),
            ],
            [
                3,
                sum([get_statistics([p0_0, p0_1, "empty"]) for p0_0, p0_1 in product(particles, particles)]),
                sum([get_statistics([p0_0, p0_1, "co"]) for p0_0, p0_1 in product(particles, particles)]),
                sum([get_statistics([p0_0, p0_1, "o"]) for p0_0, p0_1 in product(particles, particles)]),
            ]
        ]
    )

    # Format the file name.
    extension0 = ".csv"
    file_name0_ = folder_name + file_name0
    file_name1_ = cp.deepcopy(file_name0_)

    counter0 = 0
    while os.path.isfile(file_name0_ + extension0):
        file_name0_ = cp.deepcopy(file_name1_ + str(counter0))
        counter0 += 1

    # Save the file.
    with open(file_name0_ + extension0, mode="w", newline="\n") as fl0:
        writer = csv.writer(fl0, delimiter=",")
        writer.writerows(information0)


# ------------------------------------------------------------------------------
# Select Methods.
# ------------------------------------------------------------------------------


def select_cell() -> int:
    """
        Selects a cell where an action will be take place.

        :return: The numerical id of the cell where an action will take place.
    """

    # Get the total rates of the system.
    total_rate0 = get_total_rate()

    # Select a cell.
    random_rate0 = get_random_number(0.0, total_rate0)

    # Select the cell.
    cumulative0 = 0.0
    for i0, rates0 in enumerate(rates):
        cumulative0 += rates0[-1]
        if random_rate0 < cumulative0:
            return i0

    # There is an inconsistency if this point is reached.
    raise ValueError("The cell chosen to perform an action is out of range.")


def select_process(site0: int) -> int:
    """
        Selects the process to be performed.

        :param site0: The site in which the process is taking place.
    """

    # Get the total rate of the given cell.
    random_rate0 = get_random_number(0.0, rates[site0][-1])

    # Get the proper rate.
    for i0, rate0 in enumerate(rates[site0]):
        if random_rate0 < rate0:
            return i0

    # There is an inconsistency if this point is reached.
    raise ValueError("The action to be performed is out of range.")


# ------------------------------------------------------------------------------
# Update Methods.
# ------------------------------------------------------------------------------


def update_lattice(site0: int, action0: int) -> None:
    """
        Given a site in the lattice and a move to make, updates the lattice
        accordingly.

        :param site0: The numerical identifier of the site at which the action
         will take place.

        :param action0: The numerical identifier of the action that will take
         place at the site.
    """

    # Make the move:
    if action0 == 0:
        action_reaction_lh(site0)
        return

    elif action0 == 1:
        action_reaction_er(site0, "co")
        return

    elif action0 == 2:
        action_reaction_er(site0, "o")
        return

    elif action0 == 3:
        action_adsorption(site0, "co")
        return

    elif action0 == 4:
        action_desorption(site0, "co")
        return

    elif action0 == 5:
        action_diffusion(site0)
        return

    elif action0 == 6:
        action_adsorption(site0, "o")
        return

    elif action0 == 7:
        action_desorption(site0, "o")
        return

    elif action0 == 8:
        action_diffusion(site0)
        return

    raise ValueError(f"The requested action for the given site doesn't exist. Requested action identifier: {action0}.")


def update_rates() -> None:
    """
        Updates the cumulative rates of all the cells in the system.
    """

    # For all the sites.
    for i0 in range(3):
        rates[i0][0] = rate_lh(i0)  # Langmuir-Hinshelwood process.
        rates[i0][1] = rates[i0][0] + rate_er(i0, "co")  # Elay-Rideal process for carbon monoxide.
        rates[i0][2] = rates[i0][1] + rate_er(i0, "o")  # Elay-Rideal process for oxygen.
        rates[i0][3] = rates[i0][2] + rate_adsorption(i0, "co")  # Carbon monoxide adsorption.
        rates[i0][4] = rates[i0][3] + rate_desorption(i0, "co")  # Carbon monoxide desorption.
        rates[i0][5] = rates[i0][4] + rate_diffusion(i0, "co")  # Carbon monoxide diffusion.
        rates[i0][6] = rates[i0][5] + rate_adsorption(i0, "o")  # Oxygen adsorption.
        rates[i0][7] = rates[i0][6] + rate_desorption(i0, "o")  # Oxygen desorption.
        rates[i0][8] = rates[i0][7] + rate_diffusion(i0, "o")  # Oxygen diffusion.


# //////////////////////////////////////////////////////////////////////////////
# Main Program.
# //////////////////////////////////////////////////////////////////////////////


if __name__ == "__main__":
    # Run the simulation.
    run_simulation(equilibrate0=True)
    elapsed_time, steps = run_simulation(equilibrate0=False)

    print(f"Sum of probabilities: {sum(statistics)}\nDifference: {1.0 - sum(statistics)}")

    # Save to a file.
    save_statistics(elapsed_time, steps, file_name)

    # Make the plot.
    plot_pies()
