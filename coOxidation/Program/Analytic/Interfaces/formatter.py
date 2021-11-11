""" Serves as the base class to create an equation formatter."""

# ------------------------------------------------------------------------------
# Imports.
# ------------------------------------------------------------------------------

# Imports: General.
from abc import ABC, abstractmethod
from typing import Union

# ------------------------------------------------------------------------------
# Classes.
# ------------------------------------------------------------------------------


class Formatter(ABC):
    """ An abstract static class that contains base equation formatting
        functions that are implemented to generate equations in the given
        format.
    """

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Methods.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    #  Format Methods.
    # --------------------------------------------------------------------------

    @staticmethod
    @abstractmethod
    def format_constraint(constraint: tuple[tuple, tuple]) -> str:
        """ Gets the string that represents a constraint of the system in
            the given format.

            :param constraint: The variable that contains a constraint, in the
             form of an equality.

            :return: The string that represents the constraint in the given
             format.
        """
        return ""

    @staticmethod
    @abstractmethod
    def format_equation(equation: tuple, order: int = 0) -> str:
        """ Gets the string that represents an equation from a Master Equation
            in the given format.

            :param equation: The tuple that contains, in order: 1. The state for
             which the master equation will be written. 2. The dictionary of the
             states that will decay to the state for which the master equation
             will be written; where the keys are the associated decay rate
             constants for each process. Multiplicity of the states are
             included. 3. The dictionary of the states to which the state will
             decay due to the different processes; where the keys are the
             associated decay rate constants for each process. Multiplicity of
             the states are included.

            :param order: The order to which the state must be expanded. Order
             zero means the state must not be modified. Higher orders means the
             state must be mean-field expanded to the given order.

            :return: The string that represents the Master Equation in
             Mathematica format.
        """
        return ""

    @staticmethod
    @abstractmethod
    def format_final(quantities: dict) -> str:
        """ Formats the equations and the different quantities such that it
            ready to be saved in a string.

            :param quantities: A dictionary that MUST have the following format,
             with the keys: 1. "constraints": A list of strings of constraints.
             2. "equations": A list of strings of the equations. 3. "initial
             conditions": A list of the strings of the initial conditions. 4.
             "rate values": A list of the strings with the values of the
             constants. 5. "raw_states": A list with the representation of the
             raw states. For some formats the raw states may be the same as the
             regular states.

            :return: A single string of the formatted equations.
        """

        keys = ["constraints", "equations", "initial conditions", "rate values", "raw states"]
        return str(keys)

    @staticmethod
    @abstractmethod
    def format_initial_condition(state: tuple, time: Union[float, int, str] = 0.0,
                                 value: Union[float, int, str] = 0.0) -> str:
        """ Gets the string that represents a state equal to a given initial
            condition that, by default, is set to zero.

            :param state: The state whose initial condition will be.

            :param time: A parameter, that must allow a string reprsentation,
             that denotes the time of the initial condition. Set to zero by
             default.

            :param value: The value of the initial condition. Must allow a
             string representation.

            :return: The string that represents the initial condition of a
             state.
        """
        return ""

    @staticmethod
    @abstractmethod
    def format_rate(rate: str, value: Union[float, int, str] = None) -> str:
        """ Gets the string that represents a rate constant in the given format.

            :param rate: The rate to check. If in the format: "'s1'.'s2'. ...
             .'sN'..", it interprets the periods as sub-indexes of level N.

            :param value: The value to which the rate must be set.

            :return: The rate constant in the given format.
        """
        return ""

    @staticmethod
    @abstractmethod
    def format_state(state: tuple, order: int = 0, raw: bool = False):
        """ Gets the string that represents a state in the given format.

            :param state: A state in the format, ((particle0, index0), ... ,
             (particleN, indexN),).

            :param order: The order to which the state must be approximated. If
             zero, the state is NOT approximated.

            :param raw: True, if the state must be formatted as a raw state. For
             some formats, the raw state is the same as the regular state.

            :return: The state, to the given order, in the given format.
        """
        return ""

    # --------------------------------------------------------------------------
    # Get Methods.
    # --------------------------------------------------------------------------

    @staticmethod
    @abstractmethod
    def get_format_methods() -> dict:
        """ Returns a dictionary with the possible quantities to be formatted.

            :return: A dictionary with the possible quantities to be obtained by
             the formatter.
        """

        # The dictionary of the possible features to format.
        formatter_functions = {
            "constraint": Formatter.format_constraint,
            "equation": Formatter.format_equation,
            "final": Formatter.format_final,
            "initial condition": Formatter.format_initial_condition,
            "rate": Formatter.format_rate,
            "state": Formatter.format_state
        }

        return formatter_functions
