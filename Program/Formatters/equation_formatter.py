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
    # Get Formatting Functions.
    # --------------------------------------------------------------------------

    @staticmethod
    @abstractmethod
    def get_constraint(constraint):
        """ Gets the string that represents a constraint of the system in
            the given format.

            :param constraint: The variable that contains a constraint, in
            the form of an equality.

            :return constraint_string: The string that represents the constraint
            in the given format.
        """
        constraint_string = ""
        return constraint_string

    @staticmethod
    @abstractmethod
    def get_equation(equation, order=0):
        """ Gets the string that represents an equation from a Master Equation
            in the given format.

            :param equation: The variable that contains the state and its
            constituents for which to get the equation.

            :param order: The order to which the state must be expanded. Order
            zero means the state must not be modified. Higher orders means the
            state must be mean-field expanded to the given order.

            :return equation_string: The string that represents the Master
            Equation in the given format.
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
        return equation_string

    @staticmethod
    @abstractmethod
    def get_initial_condition(state, time=0, value=0):
        """ Gets the string that represents a state equal to a given initial
            condition that, by default, is set to zero.

            :param state: The state whose initial condition will be.

            :param time: A parameter, that must allow a string reprsentation,
            that denotes the time of the initial condition. Set to zero by
            default.

            :param value: The value of the initial condition. Must allow a
            string representation.

            :return constraint_string: The string that represents the initial
            condition of a state.
        """
        initial_condition_string = ""
        return initial_condition_string

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
        return rate_string

    @staticmethod
    @abstractmethod
    def get_rate_value(rate, value=0):
        """ Gets the string that represents a rate constant in the given format,
            with the given value.

            :param rate: The rate to check. If in the format,
                 "'s1'.'s2'. ... .'sN'.."
            it interprets the periods as sub-indexes of level N.

            :param value: The value of the rate. Set to zero as default.

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
        return state_string

    @staticmethod
    @abstractmethod
    def get_state_raw(state):
        """ Gets the string that represents a 'raw state' in the given format.

            :param state: A state in the format,
                ((particle0, index0), ... ,(particleN, indexN),).

            :return state_string: The state, to the given order, in the given
            format. In some cases, the raw format might be the same as the
            regular format.
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
        return state_string

    # --------------------------------------------------------------------------
    # Other Formatting Functions.
    # --------------------------------------------------------------------------

    @staticmethod
    @abstractmethod
    def join_equations(quantities):
        """ Formats the equations and the different quantities such that it
            ready to be saved in a string.

            :param quantities: A dictionary that MUST have the following
            format, with the keys AS SHOWN

                - "constraints": A list of strings of constraints.
                - "equations": A list of strings of the equations.
                - "initial conditions": A list of the strings of the initial
                  conditions.
                - "rate values": A list of the strings with the values of the
                  constants.
                - "raw_states": A list with the representation of the raw
                  states. For some formats the raw states may be the same as the
                  regular states.

            :return formatted_system: A single string of the formatted equations.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def validate_dictionary(keys0):
            """ Validates that the dictionary is consistent.
            """

            # Get ALL the keys.
            keys0_0 = [key0_1 for key0_1 in quantities.keys()]

            # Validate that it is a dictionary.
            if not isinstance(quantities, (dict,)):
                raise TypeError("The quantities parameter must be a dictionary.")

            # Validate the entries.
            if not set(keys0_0) == set(keys0):
                raise ValueError(f"The keys in the dictionary must be {keys0}."
                                 f" Current keys = {keys0_0}.")

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Dictionary keys are.
        keys = ["constraints", "equations", "initial conditions", "rate values", "raw states"]

        # Always validate first.
        validate_dictionary(keys)

        formatted_system = ""
        return formatted_system
