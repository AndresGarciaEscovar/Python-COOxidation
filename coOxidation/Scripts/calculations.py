""" Constains the functions to perform calculations. """


# ------------------------------------------------------------------------------
# Imports.
# ------------------------------------------------------------------------------

# General.
import numpy as np

# ------------------------------------------------------------------------------
# Imports.
# ------------------------------------------------------------------------------


def adsorption_rate(area0: float, pressure0: float, molar_mass0: float, temperature0: float = 300.0) -> float:
    """
        Given the parameters, returns the rate of adsorption of the given
        specie.

        :param area0: The effective surface area of adsorption, in SI units,
         i.e., m^2.

        :param pressure0: The pressure of the system in SI units, i.e., Pascals.

        :param molar_mass0: The molar mass of the specie, in SI units, i.e.,
         kg/mol.

        :param temperature0: The temperature of the environment, in SI units,
         i.e., Kelvin.

        :return The rate of adsorption of the given specie at the given TP
         value.
    """

    # Get the Boltzmann constant in Joules.
    kb0 = boltzmann_value("j")

    # The denominator of the equation.
    denominator0 = (2.0 * np.pi * molar_mass0 * kb0 * temperature0) ** 0.5

    return area0 * pressure0 / denominator0


def arrhenius_rate(energy0: float, prefactor0: float = 1e-12, temperature0: float = 300, units0: str = "j") -> float:
    """
        Given the activation energy, returns the Arrhenius rate of the system.

        :param energy0: The activation energy of the process.

        :param prefactor0: The Arrhenius rate pre-factor, in SI units, i.e,
         s^-1.

        :param temperature0: The temperature in Kelvin.

        :param units0: The units in which the Boltzmann constant must be applied.
         Currently it only allows Joules / Kelvin or eV / Kelvin.

        :return The Arrhenius rate given the parameters.
    """

    # Get the value of the Boltzmann constant.
    kb0 = boltzmann_value(units0)

    return prefactor0 * np.exp(-energy0 / (kb0 * temperature0))


def boltzmann_value(units0: str = "j"):
    """
        Returns the value of the Boltzmann constant, in the given units.

        :param units0: The units in which the Boltzmann constant must be applied.
         Currently it only allows Joules/Kelvin or electron_volts/Kelvin.

        :return The value of the Boltzmann constant in the given units.
    """

    # Normalize the string.
    tmp0 = units0.strip().lower()
    if tmp0 not in ["j", "ev"]:
        raise ValueError("The requested units for the system are not valid.")

    # Boltzmann constant.
    kb0 = 1.380649e-23 if tmp0 == "j" else 8.61733e-5

    return kb0


# ------------------------------------------------------------------------------
# Main Program.
# ------------------------------------------------------------------------------


if __name__ == "__main__":

    # --------------------------------------------------------------------------
    # For the Arrhenius rate.
    # --------------------------------------------------------------------------

    # Boltzmann constant and temperature.
    kb_eV = 8.61733e-5
    T = 600
    print("kB .T = " + f"{kb_eV*T:.6e}" + " eV")

    # --------------------------------------------------------------------------
    # For the Adsorption rate.
    # --------------------------------------------------------------------------

    # Area in meters squared.
    area = 1.347e-19
    sites = 3
    total_area = area * sites
    print("total_area = area.sites = " + f"{total_area:.6e}" + " m^2")

    # Boltzmann constant and temperature.
    kb_J = 1.380649e-23
    T = 600

    denominator = (2.0 * np.pi * kb_J * T)**0.5
    prefact = total_area / denominator

    print("pre-factor = A/(2.pi.kB.T)^1/2 = " + f"{prefact:.6e}" + " m^2/(J)^0.5")

    # --------------------------------------------------------------------------
    # Other calculations.
    # --------------------------------------------------------------------------

    print(f"{0.016 + 0.012:e}")


