""" Constains the functions to perform calculations. """


# ------------------------------------------------------------------------------
# Imports.
# ------------------------------------------------------------------------------

# General.
import copy as cp
import csv
import os
import numpy as np

from typing import Union

# User-defined
from coOxidation.Program.Utilities.context_managers import WorkingDirectorySet

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

        :return: The rate of adsorption of the given specie at the given TP
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

        :return: The Arrhenius rate given the parameters.
    """

    # Get the value of the Boltzmann constant.
    kb0 = boltzmann_value(units0)

    return prefactor0 * np.exp(-energy0 / (kb0 * temperature0))


def boltzmann_value(units0: str = "j"):
    """
        Returns the value of the Boltzmann constant, in the given units.

        :param units0: The units in which the Boltzmann constant must be applied.
         Currently it only allows Joules/Kelvin or electron_volts/Kelvin.

        :return: The value of the Boltzmann constant in the given units.
    """

    # Normalize the string.
    tmp0 = units0.strip().lower()
    if tmp0 not in ["j", "ev"]:
        raise ValueError("The requested units for the system are not valid.")

    # Boltzmann constant.
    kb0 = 1.380649e-23 if tmp0 == "j" else 8.61733e-5

    return kb0


def format_table(table0: Union[list, tuple], precision0: int = 3) -> list:
    """
        Formats the table to be saved in excel.

        :param table0: The table to be formatted.

        :param precision0: The number of significant figures with which the
         floating point numbers must be printed.

        :return: The properly formatted table to be printed.
    """

    # --------------------------------------------------------------------------
    # Auxiliary Functions.
    # --------------------------------------------------------------------------

    def column_widths(table1: list, precision1: int) -> tuple:
        """
            Gets the widths of the columns in characters.

            :param table1: The table to be formatted.

            :param precision1: The number of significant figures with which the
             floating point numbers must be printed.

            :return: The properly formatted table to be printed.
        """

        # Check that all the rows have the same number of columns.
        if not len(set(map(len, table1))) == 1:
            raise ValueError("The number of columns for the rows differ. Check before continuing.")

        # Set the maximum to the number of characters of the zeroth column.
        max1 = list(map(lambda x: len(str(x).strip()), table1[0]))

        # Do it for each row.
        for i1, row1 in enumerate(table1):
            max1[0] = max(max1[0], len(row1[0].strip()))
            max1[1:] = list(map(lambda x, y: max(x, len(f"{y:{precision1}}")), max1[1:], row1[1:]))

        return tuple(max1)

    # --------------------------------------------------------------------------
    # Implementation.
    # --------------------------------------------------------------------------

    # Convert the table to a list of lists.
    table0_ = list(map(list, table0))
    table_final0 = []

    # Get the widths.
    widths0 = column_widths(table0_, precision0)

    # Format each row.
    for i0, row0 in enumerate(table0_):
        # Reset the temporary storage.
        row0_ = list(None for _ in enumerate(row0))

        # Format the header.
        if i0 == 0:
            row0_ = [f"{column0.strip():{width0}}" for width0, column0 in zip(widths0, row0)]
            table_final0.append(cp.deepcopy(row0_))
            continue

        # Format the other entries.
        row0_[0] = f"{row0[0].strip():{widths0[0]}}"

        # Pre-format the numerical strings.
        numerical0 = list(f"{column0:.{precision0}e}" for column0 in row0[1:])
        numerical0 = list(f"{column0:{width0}}" for column0, width0 in zip(numerical0, widths0[1:]))
        row0_[1:] = numerical0

        # Append the entry.
        table_final0.append(cp.deepcopy(row0_))

    return table_final0


def get_adsorption_rates(temperature0: float = 300.0) -> list:
    """
        Gets the table with the adsorption rates.

        :param temperature0: The temperature of the system, in SI units, i.e,
         Kelvin.

         :return: The table with the adsorption rates of the system.
    """

    # Area of the substrate (m^2)
    area0 = 1.347e-19
    sites0 = 3
    total_area0 = area0 * sites0

    # Activation energies and rates.
    adsorption_header0 = [["Description", "molar mass(kg/mol)", "Pressure (PA)", "Rate in SI units."], ]

    adsorption_rates0 = []
    adsorption_rates_base0 = [
        ["CO adsorption", 2.8e-2, 0.0, 0.0],
        ["O2 adsorption", 3.2e-2, 0.0, 0.0]
    ]

    # List of pressures.
    pressures0 = [1.01 * 10.0 ** i for i in range(10)]

    # Get the range of rates for each specie.
    for i0, molar_mass0 in enumerate(adsorption_rates_base0):
        specie0 = cp.deepcopy(molar_mass0)

        # Calculate the rate for the pressure.
        for j0, pressure0 in enumerate(pressures0):
            # Get the information.
            specie0[2] = pressure0
            specie0[3] = adsorption_rate(
                area0=total_area0, pressure0=pressure0, molar_mass0=specie0[1], temperature0=temperature0
            )

            # Append the information.
            adsorption_rates0.append(cp.deepcopy(specie0))

    # Get the header.
    adsorption_header0.extend(adsorption_rates0)

    return adsorption_header0


def get_arrhenius_rates(temperature0: float = 300.0) -> list:
    """
        Gets the table with the Arrhenius rates.

        :param temperature0: The temperature of the system, in SI units, i.e,
         Kelvin.

         :return: The table with the Arrhenius rates of the system.
    """

    # The units of the system.
    units0 = "ev"

    # Arrhenius pre-factor.
    power0 = 12.5
    prefactor0 = 10 ** power0

    # Activation energies and rates.
    arrhenius_header0 = [["Description", "Activation energy(eV)", "Rate in SI units."], ]
    arrhenius_rates0 = [
        ["CO desorption", 1.10, 0.0],
        ["CO diffusion", 0.699, 0.0],
        ["O desorption", 2.87, 0.0],
        ["O diffusion", 0.608, 0.0]
    ]

    # Get the Arrhenius rates.
    for i0, energy0 in enumerate(arrhenius_rates0):
        # Get the rate.
        rate0 = arrhenius_rate(energy0=energy0[1], prefactor0=prefactor0, temperature0=temperature0, units0=units0)
        arrhenius_rates0[i0][2] = rate0

    # Add the header.
    arrhenius_header0.extend(arrhenius_rates0)

    return arrhenius_header0


def save_tables(tables0: list, filename0: str = "tables") -> None:
    """
        Save the tables in the given file. The file is saved in the "current
        working directory".

        :param tables0: The, already formatted tables, to be saved.

        :param filename0: The name of the file where the tables will be saved.
         Must be extensionless.
    """

    # Add the extension to the file.
    filename0_ = filename0 + ".csv"

    # Open the file in writing mode.
    with open(filename0_, newline="\n", mode="w") as fl0:
        writer = csv.writer(fl0, delimiter=",")

        # Save all the tables.
        for table0 in tables0:
            writer.writerows(table0)

# ------------------------------------------------------------------------------
# Main Program.
# ------------------------------------------------------------------------------


if __name__ == "__main__":

    # --------------------------------------------------------------------------
    # System general parameters.
    # --------------------------------------------------------------------------

    # The temperature of the system (K).
    temperature = 600

    # --------------------------------------------------------------------------
    # Get the rates.
    # --------------------------------------------------------------------------

    # Get the Arrhenius rates.
    arrhenius_rates = get_arrhenius_rates(temperature0=600.0)

    # Get the Adsorption rates.
    adsorption_rates = get_adsorption_rates(temperature0=600.0)

    # --------------------------------------------------------------------------
    # Format the Rates.
    # --------------------------------------------------------------------------

    # Format the tables.
    adsorption_rates = format_table(adsorption_rates, precision0=6)
    arrhenius_rates = format_table(arrhenius_rates, precision0=6)

    # --------------------------------------------------------------------------
    # Save the Information.
    # --------------------------------------------------------------------------

    # Save the tables.
    directory = os.path.dirname(__file__)
    with WorkingDirectorySet(os.path.dirname(__file__)) as fl:
        save_tables([arrhenius_rates, adsorption_rates], "tables")
