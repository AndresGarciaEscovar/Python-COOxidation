""" Stores the references to ALL the formatters, so they can be easily accessed.
"""

# Imports: User-defined imports.
from .formatter_latex import LaTeXFormatter
from .formatter_mathematica import MathematicaFormatter


class GetFormatter:
    """ A static class that gets the requested formatter.
    """

    # Dictionary of formatters.
    _FORMATTERS = {
        "latex": LaTeXFormatter,
        "mathematica": MathematicaFormatter
    }

    @staticmethod
    def get_formatter(requested_formatter):
        """ From the allowed list of formatters, gets a formatter.

            :param requested_formatter: Gets the requested formatter, provided
            that it is in the list.

            :return formatter: The requested formatter.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def validate_formatter():
            """ Checks that the requested formatter exists in the list possible
                values. The formatter is NOT case sensitive, e.g., 'A' = 'a'.
            """

            # Format the formatter properly.
            requested_formatter0 = requested_formatter.strip().lower()

            # Get the list of possible formatters.
            keys0 = GetFormatter._FORMATTERS.keys()

            # Check that the requested formatter is in the dictionary.
            if requested_formatter0.strip().lower() not in keys0:
                raise ValueError(f"The requested formatter must be in the list: {keys0}."
                                 f" Requested formatter: {requested_formatter0}.")

            return requested_formatter0

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Validate the formatter.
        requested_formatter = validate_formatter()

        return GetFormatter._FORMATTERS[requested_formatter]

    @staticmethod
    def get_formatter_list():
        """ Returns a list with the strings of the available formatters.

            :return formatters: A list with the available formatters.
        """

        # Get the list.
        formatters = [key for key in GetFormatter._FORMATTERS.keys()]

        return formatters
