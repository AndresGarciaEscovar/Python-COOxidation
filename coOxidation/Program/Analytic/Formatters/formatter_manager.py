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
    # Public Interface.
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

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def validate_formatter(formatter0: str):
            """ Checks that the requested formatter exists in the list possible
                values. The formatter is NOT case sensitive, e.g., 'A' = 'a'.

                :param formatter0: The string that represents the requested
                 formatter static class.
            """

            # Format the formatter properly.
            formatter0_ = formatter0.strip().lower()

            # Get the list of possible formatters.
            keys0 = FormatterManager._FORMATTERS.keys()

            # Check that the requested formatter is in the dictionary.
            if formatter0_ not in keys0:
                raise ValueError(f"The requested formatter must be in the list: {keys0}."
                                 f" Requested formatter: {formatter0_}.")

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        # Validate the formatter.
        validate_formatter(formatter)
        
        # Delete blank spaces and make it lower case.
        formatter_ = formatter.strip().lower()

        return FormatterManager._FORMATTERS[formatter_]

    @staticmethod
    def get_formatter_list() -> list:
        """ Returns a list with the strings of the available formatters.

            :return: A list with the available formatters.
        """
        return [key for key in FormatterManager._FORMATTERS.keys()]
