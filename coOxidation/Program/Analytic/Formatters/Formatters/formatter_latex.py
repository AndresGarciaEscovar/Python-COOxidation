""" Equation formatter for LaTeX."""

# ------------------------------------------------------------------------------
# Imports.
# ------------------------------------------------------------------------------

# Imports: General.
from typing import Union

# Imports: User-defined.
from coOxidation.Program.Analytic.Interfaces.formatter import Formatter

# ------------------------------------------------------------------------------
# Classes.
# ------------------------------------------------------------------------------


class LaTeXFormatter(Formatter):
    """ A static class that contains equation formatting functions that are
        generated by the equation generator.
    """

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Public Interface.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Get Methods.
    # --------------------------------------------------------------------------

    @staticmethod
    def get_constraint(constraint: tuple) -> str:
        """ Gets the string that represents a constraint of the system in
            LaTeX format.

            :param constraint: The variable that contains the constraint, in the
             form of equalities. It must be a tuple with two entries such that
             the first entry represents the right-hand side of the equation,
             i.e, a single state, and the second entry the left-hand side of the
             equation, i.e., a list of a single, or multiple, states.

            :return: The string that represents the constraint in LaTeX format.
        """

        # Get the lowest order state.
        low_state = LaTeXFormatter.get_state(constraint[0])

        # Get the other states.
        other_states = list(map(LaTeXFormatter.get_state, constraint[1]))

        # Join the states.
        other_states = " + ".join(other_states)

        # Join the strings.
        constraint_ = low_state + " = " + other_states

        return constraint_

    @staticmethod
    def get_equation(equation: tuple, order: int = 0) -> str:
        """ Gets the string that represents an equation from a Master Equation
            in LaTeX format.

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

            :return: The string that represents the Master Equation in LaTeX
             format.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def format_create_decay(key0, create_states0, decay_states0):
            """ Given the decay states and the key, it formats the string of
                decay states.

                :param key0: The key that is being formatted.

                :param create_states0: The create states associated with the
                key.

                :param decay_states0: The decay states associated with the key.

                :return create_decay_string0: The string that represents the
                specific term in the equation.
            """

            # Get the string representation of the key.
            string_key = LaTeXFormatter.get_rate(key0)

            # Join the states in the create states list.
            create_string0 = "+".join(create_states0)

            # Join the states in the decay states list.
            decay_string0 = "-" + "-".join(decay_states0)

            # Join the strings.
            create_decay_string0 = f"+{string_key} \\left(" + create_string0 + decay_string0 + f"\\right)"

            return create_decay_string0

        def format_create_decay_single(key0, states0, decay=False):
            """ Given the decay states and the key, it formats the string of
                decay states.

                :param key0: The key that is being formatted.

                :param states0: The create/decay states associated with the key.

                :param decay: True, if the requested states to be added are
                decay states. False, otherwise, i.e., create states.

                :return create_decay_string: The string that represents the
                specific term in the equation.
            """

            # ------------------------------------------------------------------
            # Auxiliary functions.
            # ------------------------------------------------------------------

            def get_prefactor(state_string1):
                """ Given a state string, it returns the string of the numerical
                    coefficient, and the stripped state.
                """

                # Auxiliary variables.
                j = 0
                tmp_string1 = ""

                # Every character in the string.
                for character in state_string1:
                    # If the character is not a number.
                    if not str.isnumeric(character):
                        break

                    # Add the character to the string.
                    tmp_string1 += character

                    # Add one to the counter.
                    j += 1

                # Format the string properly.
                state_string1 = state_string1[j:] if j < len(state_string1) else state_string1

                # THIS SHOULD NOT HAPPEN!
                return tmp_string1, state_string1

            # ------------------------------------------------------------------
            # Implementation.
            # ------------------------------------------------------------------

            # Get the string representation of the key.
            string_key = LaTeXFormatter.get_rate(key0)

            # Initialize the string and the negative sign as needed.
            create_decay_string0 = "-" if decay else "+"

            # If there is only one state.
            if len(states0) == 1:
                # Get the prefactor and state.
                prefactor0, create_decay_string0_0 = get_prefactor(states0[0])

                # Join the string.
                create_decay_string0 += f"{prefactor0} {string_key} {create_decay_string0_0}"

            else:
                # Join the states.
                create_decay_string0 += f"{string_key} \\left(" + "+".join(states0) + "\\right)"

            return create_decay_string0

        def format_equation_string(equation_string0):
            """ Formats the equation string further to include spaces for
                readability.

                :param equation_string0: The string to be formatted.

                :return equation_string0_0: The properly formatted string.
            """

            # Strip all the leading and trailing spaces.
            equation_string0_0 = equation_string0.strip()

            # Save the negative character if needed.
            first_character0 = "-" if equation_string0[0] == "-" else ""

            # Determine if there is a positive or negative sign at the start.
            delete_first0 = equation_string0_0[0] == "-" or equation_string0_0[0] == "+"
            equation_string0_0 = equation_string0_0[1:] if delete_first0 else equation_string0_0

            # Strip all the leading and trailing spaces, again.
            equation_string0_0 = equation_string0_0.strip()

            # Space the positive and negative signs correctly.
            equation_string0_0 = " + ".join(equation_string0_0.split("+"))
            equation_string0_0 = " - ".join(equation_string0_0.split("-"))

            # Add the first character.
            equation_string0_0 = first_character0 + equation_string0_0

            return equation_string0_0

        def format_state_multiplicity(state0):
            """ Returns the state string, properly formatted, multiplied by its
                multiplicity.

                :param state0: A 2-tuple of the state with its multiplicity.

                :return: The state string, properly formatted, multiplied by its
                multiplicity
            """

            # Check that the state is an iterable of length 2.
            if not len(state0) == 2:
                raise ValueError("To properly format the state it must be a"
                                 " tuple of lenght 2.")

            # Get the multiplicity.
            state0_0 = str(state0[1]) if state0[1] > 1 else ""

            # Get the state representation.
            state0_0 += LaTeXFormatter.get_state(state0[0], order)

            return state0_0

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
            keys0_1 = set(key0 for key0 in equation[1].keys())
            keys0_2 = set(key0 for key0 in equation[2].keys())

            # If their keys are different.
            if not keys0_1 == keys0_2:
                raise ValueError(f"Both Dictionaries must have the same keys"
                                 f" Keys for equation[1]: {keys0_1},"
                                 f" Keys for equation[2]: {keys0_2}.")

        def validate_state(state0):
            """ Validates that the state is given in the proper format.

                :param state0: A state that must be in the format,
                    ((particle0, index0), ... ,(particleN, indexN),).
            """

            # Check that it is a tuple.
            if not isinstance(state0, (tuple,)):
                raise TypeError(f"The state parameter must be a tuple. Current type: {type(state0)}")

            # Check that the elements are tuples.
            for j0, substate0 in enumerate(state0):
                # Check that it is a tuple.
                if not isinstance(substate0, (tuple,)):
                    raise TypeError(f"The substates of a state must be tuple."
                                    f" State = {state0}, Substate Entry = {j0},"
                                    f" Substate = {substate0}, Current type: {type(substate0)}"
                                    )

                # Of length 2.
                elif not len(substate0) == 2:
                    raise TypeError(f"The substates of a state must be tuple of length 2. "
                                    f" State = {state0}, Substate Entry = {j0},"
                                    f" Substate = {substate0},  Current length of substate: {len(substate0)}."
                                    )

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Validate the equation.
        validate_equation()

        # Auxiliary variables.
        keys = tuple(key for key in equation[1].keys())

        # Get the differential form.
        diff_state = "\\frac{d" + LaTeXFormatter.get_state(equation[0]) + "}{dt} ="

        # The string where the equation will be stored.
        equation_string = ""

        # For every key.
        for key in keys:
            # Get the create states representations.
            create_states = [format_state_multiplicity(state) for state in equation[1][key]]

            # Get the decay states representations.
            decay_states = [format_state_multiplicity(state) for state in equation[2][key]]

            # If there are decay states but no creation states.
            if len(decay_states) > 0 and len(create_states) == 0:
                # Get the decay string.
                decay_string = format_create_decay_single(key, decay_states, decay=True)

                # Add to the equation string.
                equation_string += decay_string

            # If there are no decay states, but there are creation states.
            elif len(decay_states) == 0 and len(create_states) > 0:
                # Get the create string.
                create_string = format_create_decay_single(key, create_states, decay=False)

                # Add to the equation string.
                equation_string += create_string

            # If there are both decay states and creation states.
            elif len(decay_states) > 0 and len(create_states) > 0:
                # Get the decay and create string.
                create_decay_string = format_create_decay(key, create_states, decay_states)

                # Add to the equation string.
                equation_string += create_decay_string

        # Format the string properly.
        equation_string = format_equation_string(equation_string)

        # Join the strings.
        equation_string = diff_state + equation_string

        return equation_string

    @staticmethod
    def get_initial_condition(state: tuple, time: Union[float, str] = 0.0, value: Union[float, str] = 0.0) -> str:
        """ Gets the string that represents a state equal to a given initial
            condition that, by default, is set to zero.

            :param state: The state whose initial, at the given time, takes a
             given value.

            :param time: A floating point number, or a string, that represents
             the initial time at which the initial condition is set. Set to zero
             by default.

            :param value: A floating point number, or a string, that represents
             the initial condition at the initial time. Set to zero by default.

            :return: The string that represents the initial condition of a
             time-dependent state.
        """
        # Add a sub-index to the state.
        initial_condition = LaTeXFormatter.get_state(state) + "_{t_{0} =" + f"{str(time)}" + "} = " + f"{value}"

        return initial_condition

    @staticmethod
    def get_rate(rate: str) -> str:
        """ Gets the string that represents a rate constant in LaTeX format.

            :param rate: The rate string to turn into a LaTeX variable. If it is
             in the format: "'s1'.'s2'. ... .'sN'..", it interprets the periods
             as sub-indexes of level N.

            :return: The rate constant in LaTeX format.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def validate_rate(rate0: str) -> None:
            """ Validates that the rate is given in the proper format.

                :param rate0: The representation of the rate that must be a
                 string.
            """

            # Check that it is a string.
            if not isinstance(rate0, (str,)):
                raise TypeError(f"The rate parameter must be string. Current type: {type(rate0)}")

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        # Validate the form of the rate.
        validate_rate(rate)

        # Split the string.
        rate_string = rate.split(".")

        # Join the string properly.
        rate_string = "_{".join(rate_string) + "}" * (len(rate_string) - 1)

        return rate_string

    @staticmethod
    def get_rate_value(rate: str, value: float = 0.0) -> str:
        """ Gets the string that represents a rate constant in LaTeX format,
            with the given value.

            :param rate: The rate string to turn into a LaTeX variable. If it is
             in the format: "'s1'.'s2'. ... .'sN'..", it interprets the periods
             as sub-indexes of level N.

            :param value: The numerical value of the rate. Set to zero as
             default.

            :return: The rate constant in LaTeX format.
        """

        # Get the rate representation.
        rate_string = LaTeXFormatter.get_rate(rate)

        # Set the value.
        rate_string += f" = {value}"

        return rate_string

    @staticmethod
    def get_state(state: tuple, order: int = 0) -> str:
        """ Gets the string that represents a state in LaTeX format.

            :param state: A tuple that represents the state in the format,
             ((particle0, index0), ... ,(particleN, indexN),).

            :param order: The order to which the state must be approximated.
             If zero, the state is NOT approximated.

            :return: The state, to the given order, in LaTeX format.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def format_component(component0: tuple) -> str:
            """ Formats an entry of a state properly.

                :param component0: The component of a state to be formatted.

                :return: The component of a state in string format.
            """

            # Get the representation of the single component.
            state0 = f"{component0[0]}" + "_{" + f"{component0[1]}" + "}"

            return state0

        def get_denominator(numerator0: list) -> list:
            """ Gets the denominator for the equations. This is obtained from
                the numerator.

                :param numerator0: A list of the states that represent the
                 numerator of the system.

                :return: The list of states that is generated from the
                 numerator, that will go in the denominator.
            """

            # Auxiliary variables.
            denominator0 = []

            # For every ith state.
            for i0, state0_0 in enumerate(numerator0):
                # For every jth state.
                for j0, state0_1 in enumerate(numerator0):
                    # The denominator will be the intersection of the states.
                    if j0 <= i0:
                        continue

                    # Get the intersecting states.
                    intersection0 = list(set(state0_0).intersection(set(state0_1)))

                    # If there are intersecting states.
                    if len(intersection0) > 0:
                        # Extend the list if there is intersection.
                        denominator0.append(intersection0)

            return denominator0

        def get_numerator(state0: tuple, order0: int) -> list:
            """ Gets the split state to the nth order, using a mean-field
                approximation.

                :param state0: A tuple of tuples that represents the state.

                :param order0: The order to which the equation must be
                 approximated.

                :return: A list of the states that is generated from an nth
                 order mean field approximation.
            """

            # The list where the numerator terms will be stored.
            numerator0 = []

            # For every index in the state.
            for i0, _ in enumerate(state0):
                # Cannot generate more states.
                if i0 + order0 > len(state0):
                    break

                # Append the substate.
                numerator0.append(state0[i0: i0 + order0])

            return numerator0

        def get_state(state0: tuple) -> str:
            """ Given a state it returns the string representation.

                :param state0: A tuple of 2-tuples.

                :return: The CLOSED reprensentation of a state.
            """

            # Get list of the formatted string for each sub-state.
            substate0 = list(map(format_component, state0))

            # Get the exact representation of the state.
            state0_ = "\\left<" + ",".join(substate0) + "\\right>"

            return state0_

        def validate_state(state0: tuple, order0: int) -> None:
            """ Validates that the state is given in the proper format.

                :param state0: The state to be formatted.

                :param order0: The order to which the state must be
                 approximated.
            """

            # The order must a number greater than zero.
            if not 0 <= order0:
                raise ValueError(f"The order must be a positive value. Current Value {order0}.")

            # Check that it is a tuple.
            if not isinstance(state0, (tuple,)):
                raise TypeError(f"The state parameter must be a tuple. Current type: {type(state0)}")

            # Check that the elements are tuples.
            for i0, substate0 in enumerate(state0):
                # Check that it is a tuple.
                if not isinstance(substate0, (tuple,)):
                    raise TypeError(f"The substates of a state must be tuple."
                                    f" State = {state0}, Substate Entry = {i0},"
                                    f" Substate = {substate0}, Current type: {type(substate0)}"
                                    )

                # Of length 2.
                elif not len(substate0) == 2:
                    raise TypeError(f"The substates of a state must be tuple of length 2. "
                                    f" State = {state0}, Substate Entry = {i0},"
                                    f" Substate = {substate0},  Current length of substate: {len(substate0)}."
                                    )

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Always validate the state.
        validate_state(state, order)

        # If the requested order is zero.
        if order == 0 or order >= len(state):
            # Get the string representation of the state.
            state_string = get_state(state)

            return state_string

        # Split the state in the requested order.
        numerator = get_numerator(state, order)

        # Get the denominator states.
        denominator = get_denominator(numerator)

        # ----------------------------------------------------------------------
        # Get the strings for the numerator and denominator.
        # ----------------------------------------------------------------------

        # Strings for the numerator states.
        numerator = list(map(get_state, numerator))

        # Format the string.
        state_string = "".join(numerator)

        # Get the denominator strings.
        if len(denominator) > 0:
            # Strings for the denominator states.
            denominator = list(map(get_state, denominator))

            # Format the string.
            state_string = "\\frac{" + state_string + "}{" + "".join(denominator) + "}"

        return state_string

    @staticmethod
    def get_state_raw(state: tuple) -> str:
        """ Gets the string that represents a 'raw state' in LaTeX format.

            :param state: A state in the format ((particle0, index0), ... ,
             (particleN, indexN),).

            :return: The string that represents the state, to the given order,
             in LaTeX format. In some cases, the raw format might be the
             same as the regular format.
        """

        # Get the string representation of the state.
        state_string = LaTeXFormatter.get_state(state)

        return state_string

    # --------------------------------------------------------------------------
    # Join Methods.
    # --------------------------------------------------------------------------

    @staticmethod
    def join_equations(quantities: dict) -> str:
        """ Formats the equations and the different quantities such that they
            are ready to be saved in a string.

            :param quantities: A dictionary that MUST have the following
             keys and variables: 1. "constraints": A list of strings of
             constraints. 2. "equations": A list of strings of the equations.
             3. "initial conditions": A list of the strings of the initial
             conditions. 4. "rate values": A list of the strings with the values
             of the constants. 5. "raw_states": A list with the representation
             of the raw states. For some formats the raw states may be the same
             as the regular states.

            :return: A string of the formatted equations.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def format_constraints(quantities0: dict) -> str:
            """ Formats the constraints of the system such that they are printed
                as a collection of LaTeX functions.

                :param quantities0: The dictionary with the quantities to be
                 formatted.

                :return: A string that contains the collection of constraints,
                 formatted for LaTeX.
            """

            # Join the constrains list.
            constraint_list0 = "\n".join(quantities0["constraints"])

            return constraint_list0

        def format_equations(quantities0: dict) -> str:
            """ Formats the equations such that they are printed as a
                LaTeX list.

                :param quantities0: The dictionary with the quantities to be
                 formatted.

                :return: The string that represents a list of LaTeX differential
                 equations, with the initial conditions.
            """

            # Join the list of equations with the initial conditions.
            equations0 = "\n".join(quantities0["equations"]) + "\n" + "\n".join(quantities0["initial conditions"])

            # Format spacing of the equality sign properly.
            equations0 = "= ".join(equations0.split("="))

            # Format it so that it reads easily when exported to LaTeX.
            equations0 = equations0[:-1] if equations0[-1] == "\n" else equations0

            return equations0

        def format_rates(quantities0: dict) -> str:
            """ Formats the rates such that they are printed as a collection of
                LaTeX equalities.

                :param quantities0: The dictionary with the quantities to be
                 formatted.

                :return: The list of formatted rates, ready to be placed in the
                 system.
            """

            # Join the list of rates.
            rates_list0 = "\n".join(quantities0["rate values"])

            return rates_list0

        def format_raw_states(quantities0: dict) -> str:
            """ Formats the raw states, i.e., the names of the states without
                any time dependence, printed as a LaTeX list.

                :param quantities0: The dictionary with the quantities to be
                 formatted.

                :return: The raw states printed as a LaTeX list.
            """

            # Join the list of equations with the initial conditions.
            raw_states0 = "\n".join(quantities0["raw states"])

            return raw_states0

        def validate_dictionary(quantities0: dict, keys0: list) -> None:
            """ Validates that the dictionary is consistent.

                :param quantities0: The dictionary with the quantities to be
                 formatted.

                :param keys0: The keys in the dictionary.
            """

            # Get ALL the keys.
            keys0_ = [key0_1 for key0_1 in quantities0.keys()]

            # Validate that it is a dictionary.
            if not isinstance(quantities0, (dict,)):
                raise TypeError("The quantities parameter must be a dictionary.")

            # Validate the entries.
            if not set(keys0_) == set(keys0):
                raise ValueError(f"The keys in the dictionary must be {keys0}."
                                 f" Current keys = {keys0_}.")

        # //////////////////////////////////////////////////////////////////////
        # Implementation
        # //////////////////////////////////////////////////////////////////////

        # ----------------------------------------------------------------------
        # Get the keys and validate the dictionary.
        # ----------------------------------------------------------------------

        # Dictionary keys are.
        keys = ["constraints", "equations", "initial conditions", "rate values", "raw states"]

        # Always validate first.
        validate_dictionary(quantities, keys)

        # ----------------------------------------------------------------------
        # Get strings of the different quantities.
        # ----------------------------------------------------------------------

        # Get the constraints string.
        constraints_list = format_constraints(quantities)

        # Get the formatted equations and initial conditions string.
        equations_list = format_equations(quantities)

        # Get the rates list string.
        rates_list = format_rates(quantities)

        # Get the list of raw variables string.
        raw_states_list = format_raw_states(quantities)

        # ----------------------------------------------------------------------
        # Join the final string.
        # ----------------------------------------------------------------------

        # Join the strings in the proper order string.
        formatted_system = rates_list + "\n" + raw_states_list + "\n" + equations_list + "\n" + constraints_list

        return formatted_system
