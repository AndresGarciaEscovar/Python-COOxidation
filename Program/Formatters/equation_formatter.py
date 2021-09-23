""" Serves as the base class to create an equation formatter.
"""

# Imports: General.
from abc import ABC, abstractmethod


class EquationFormatter(ABC):
    """ An abstract static class that contains base equation formatting
        functions that are implemented to generate equations in the given
        format.
    """

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Public Interface.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Get All Functions.
    # --------------------------------------------------------------------------

    @staticmethod
    def get_format_functions():
        """ Returns a dictionary with the possible quantities to be formatted.

            :return formatter_functions: A dictionary with the possible
            quantities to be obtained by the formatter.
        """

        # The dictionary of the possible features to format.
        formatter_functions = {
            "equation": EquationFormatter.get_equation,
            "rate": EquationFormatter.get_rate,
            "state": EquationFormatter.get_state
        }

        return formatter_functions

    # --------------------------------------------------------------------------
    # Formatting Functions.
    # --------------------------------------------------------------------------

    @staticmethod
    @abstractmethod
    def get_equation(equation, order=0):
        """ Gets the string that represents a state in the given format.

            :param equation: The variable that contains the state and its
            constituents for which to get the equation.

            :param order: The order to which the state must be expanded. Order
            zero means the state must not be modified. Higher orders means the
            state must be mean-field expanded to the given order.

            :return equation_string: The string that represents the state in the
            given format.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def validate_equation():
            """ Validates that the equation is given in the proper format.
                    (state, create states dictonary, decay states dictionary)
            """

            # Check that the equation is tuple.
            if not isinstance(equation, (tuple,)):
                raise TypeError(f"The equation must be a tuple. Current type: {type(equation)}.")

            # Of length 3.
            elif not len(equation) == 3:
                raise TypeError(f"The equation must be a tuple of three entries."
                                f" Current type: {type(equation)}, Length = {len(equation)}"
                                )

            # Validate that the zeroth entry is a state.
            validate_state(equation[0])

            # Validate that the first and second entries are dictionaries.
            if not (isinstance(equation[1], (dict,)) and isinstance(equation[2], (dict,))):
                raise TypeError("The two last entries of the tuple must be dictionaries. "
                                f" Dictionary entry [1] = {type(equation[1])},"
                                f" Dictionary entry [2] = {type(equation[2])}."
                                )

            # Get the keys to the dictionaries.
            keys1 = set(key for key in equation[1].keys())
            keys2 = set(key for key in equation[2].keys())

            # If their keys are different.
            if not keys1 == keys2:
                raise ValueError(f"Both Dictionaries must have the same keys"
                                 f" Keys for equation[1]: {keys1},"
                                 f" Keys for equation[2]: {keys2}.")

        def validate_state(state):
            """ Validates that the state is given in the proper format.

                :param state: A state that must be in the format,
                    ((particle0, index0), ... ,(particleN, indexN),).
            """

            # Check that it is a tuple.
            if not isinstance(state, (tuple,)):
                raise TypeError(f"The state parameter must be a tuple. Current type: {type(state)}")

            # Check that the elements are tuples.
            for j, substate in enumerate(state):
                # Check that it is a tuple.
                if not isinstance(substate, (tuple,)):
                    raise TypeError(f"The substates of a state must be tuple."
                                    f" State = {state}, Substate Entry = {j},"
                                    f" Substate = {substate}, Current type: {type(substate)}"
                                    )

                # Of length 2.
                elif not len(substate) == 2:
                    raise TypeError(f"The substates of a state must be tuple of length 2. "
                                    f" State = {state}, Substate Entry = {j},"
                                    f" Substate = {substate},  Current length of substate: {len(substate)}."
                                    )

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Always validate the equation.
        validate_equation()

        equation_string = ""
        pass
        return equation_string

    @staticmethod
    @abstractmethod
    def get_rate(rate):
        """ Gets the string that represents a rate constant in the given format.

            :param rate: The rate to check. If in the format,
                 "'s1'.'s2'. ... .'sN'.."
            it interprets the periods as sub-indexes of level N.

            :return rate_string: The rate constant in the given format.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def validate_rate():
            """ Validates that the rate is given in the proper format.

                :param rate: The rate to check. If in the format,
                     "'s1'.'s2'. ... .'sN'.."
                it interprets the periods as sub-indexes of level N.
            """

            # Check that it is a string.
            if not isinstance(rate, (str,)):
                raise TypeError(f"The rate parameter must be string. Current type: {type(rate)}")

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Always validate the equation.
        validate_rate()

        rate_string = ""
        pass
        return rate_string

    @staticmethod
    @abstractmethod
    def get_state(state, order=0):
        """ Gets the string that represents a state in the given format.

            :param state: A state in the format,
                ((particle0, index0), ... ,(particleN, indexN),).

            :param order: The order to which the state must be approximated.
            If zero, the state is NOT approximated.

            :return state_string: The state, to the given order, in LaTeX
            format.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def validate_state():
            """ Validates that the state is given in the proper format.

                :param state: A state that must be in the format,
                    ((particle0, index0), ... ,(particleN, indexN),).
            """

            # Check that it is a tuple.
            if not isinstance(state, (tuple,)):
                raise TypeError(f"The state parameter must be a tuple. Current type: {type(state)}")

            # Check that the elements are tuples.
            for j, substate in enumerate(state):
                # Check that it is a tuple.
                if not isinstance(substate, (tuple,)):
                    raise TypeError(f"The substates of a state must be tuple."
                                    f" State = {state}, Substate Entry = {j},"
                                    f" Substate = {substate}, Current type: {type(substate)}"
                                    )

                # Of length 2.
                elif not len(substate) == 2:
                    raise TypeError(f"The substates of a state must be tuple of length 2. "
                                    f" State = {state}, Substate Entry = {j},"
                                    f" Substate = {substate},  Current length of substate: {len(substate)}."
                                    )

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Always validate the state.
        validate_state()

        state_string = ""
        pass
        return state_string
