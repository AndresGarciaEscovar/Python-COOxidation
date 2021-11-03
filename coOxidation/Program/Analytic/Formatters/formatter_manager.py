""" Stores the references to ALL the formatters, so they can be easily accessed.
"""

# ------------------------------------------------------------------------------
# Imports.
# ------------------------------------------------------------------------------

# Imports: General.
from typing import Union

# Imports: User-defined imports.
from coOxidation.Program.Analytic.Formatters.Formatters.formatter_latex import LaTeXFormatter
from coOxidation.Program.Analytic.Formatters.Formatters.formatter_mathematica import MathematicaFormatter

# ------------------------------------------------------------------------------
# Classes.
# ------------------------------------------------------------------------------


class FormatterManager:
    """ A static class that gets the requested formatter."""

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Variables and Constants.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Constants.
    # --------------------------------------------------------------------------

    # Dictionary of formatters.
    _FORMATTERS = {
        "latex": LaTeXFormatter,
        "mathematica": MathematicaFormatter
    }

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Methods.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Get Methods.
    # --------------------------------------------------------------------------

    @staticmethod
    def get_formatter(formatter: str) -> Union[LaTeXFormatter, MathematicaFormatter]:
        """ From the allowed list of formatters, gets a formatter.

            :param formatter: The string that represents the requested formatter
             static class.

            :return:  The requested formatter.
        """

        formatter_ = formatter.strip().lower()
        return FormatterManager._FORMATTERS[formatter_]

    @staticmethod
    def get_formatter_list() -> list:
        """ Returns a list with the strings of the available formatters.

            :return: A list with the available formatters.
        """
        return [key for key in FormatterManager._FORMATTERS.keys()]
