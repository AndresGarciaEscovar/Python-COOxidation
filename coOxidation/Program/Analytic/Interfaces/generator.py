""" The class that serves as the base to create an equation generator.
"""

# ------------------------------------------------------------------------------
# Imports.
# ------------------------------------------------------------------------------

# Imports: General.
import copy as cp
import itertools
import numpy as np

from abc import ABC, abstractmethod
from itertools import product
from typing import Iterable, Union

# ------------------------------------------------------------------------------
# Classes.
# ------------------------------------------------------------------------------


class Generator(ABC):
    """ Generates the differential equations in different formats.

        - self.constraints: The list where the constraint equations will be
          saved.

        - self.equations: The list where the equations will be saved.

        - self.sites_number: The number of sites the system has.

        - self.states: the unique states a site can take at a given time.
    """

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Constants and Variables.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Getters, Setters and Deleters.
    # --------------------------------------------------------------------------

    @property
    def constraints(self) -> list:
        """ Gets the equations of the system.

            :return: The list with the constraints of the system.
        """
        return self.__constraints

    @constraints.setter
    def constraints(self, constraints: list) -> None:
        """ Sets the constraints of the system.

            :param constraints: The constraints of the system.
        """
        self.__constraints = list(constraints)

    @constraints.deleter
    def constraints(self) -> None:
        """ Cannot delete this variable."""
        raise AttributeError("Cannot delete the constraints.")

    # --------------------------------------------------------------------------

    @property
    def equations(self) -> list:
        """ Gets the equations of the system.

            :return: A list with the equations of the system.
        """
        return self.__equations

    @equations.setter
    def equations(self, equations: list) -> None:
        """ Sets the list of equations of the system.

            :param equations: A list with the equations system.
        """
        self.__equations = list(equations)

    @equations.deleter
    def equations(self) -> None:
        """ Cannot delete this variable."""
        raise AttributeError("Cannot delete the equations.")

    # --------------------------------------------------------------------------

    @property
    def sites_number(self) -> int:
        """ Gets the number of sites in the lattice.

            :return: Returns the number of sites that the system has.
        """
        return self.__sites_number

    @sites_number.setter
    def sites_number(self, sites_number: int) -> None:
        """ Sets the number of sites in the lattice.

            :param sites_number: Sets the number of sites that the system has.
        """

        # Set the sites number.
        self.__sites_number = int(sites_number)

        # Validate the sites number.
        self._validate_sites_number()

    @sites_number.deleter
    def sites_number(self) -> None:
        """ Cannot delete this variable."""
        raise AttributeError("Cannot delete the sites_number.")

    # --------------------------------------------------------------------------

    @property
    def states(self) -> tuple:
        """ Gets the possible states a system can have.
        """
        return cp.deepcopy(self.__states)

    @states.setter
    def states(self, states: Iterable) -> None:
        """ Sets the possible states a system can have.

            :param states: An iterable that contains the unique states a site
             can take at a given time.
        """
        self.__states = tuple(map(str, set(states)))

    @states.deleter
    def states(self) -> None:
        """ Cannot delete this variable."""
        raise AttributeError("Cannot delete the states.")

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Public Interface.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Get Methods.
    # --------------------------------------------------------------------------

    def get_nth_order_equations(self, order: int = 0, display: bool = False) -> None:
        """ Gets the nth order equations for the system in the format: (state,
            decay_states, create_states); where "decay_states" and
            "create_states" are the dictionaries that contain the states to
            which the "state" decays, or the states from where the state is
            created, due to the different process in the system.

            Saves the equations to the variables for them to be processed later.

            :param order: The lowest order to which the equations are to be
             calculated.

            :param display: If the table of equations must be displayed in the
             console.
        """

        # ----------------------------------------------------------------------
        # Empty the equations list.
        # ----------------------------------------------------------------------

        # Empty the equations list.
        self.equations = []

        # ----------------------------------------------------------------------
        # Get the states to perform the calculation.
        # ----------------------------------------------------------------------

        # Get the lowest order states.
        states_left = self._get_states_left(order)

        # Get the other involved states.
        states_right = self._get_states_right(order)

        # ----------------------------------------------------------------------
        # Get the decay states.
        # ----------------------------------------------------------------------

        # An alias to the function to make it shorter.
        function = self._get_decay_states

        # Dictionary of processes.
        processes = self._get_process_functions()

        # Get the decay states for each of the right-hand states.
        decay_states = list(map(lambda x: function(x, processes), states_right))

        # ----------------------------------------------------------------------
        # Get the create and decay state(s) of the left-hand states.
        # ----------------------------------------------------------------------

        # Get both the decay and creation dictionaries for each state.
        for state_left in states_left:
            # Get the list of create states.
            dictionary_create = self._get_products_create(state_left, decay_states)

            # Get the list of decay states.
            dictionary_decay = self._get_products_decay(state_left, decay_states)

            # Append it to the equations.
            self.equations.append((state_left, dictionary_create, dictionary_decay,))

        # Get the constraints.
        self._get_constraints(states_left)

        # If the equations must be displayed.
        if display:
            # Print the equations.
            self.print_equation_states()

    # --------------------------------------------------------------------------
    # Print Methods.
    # --------------------------------------------------------------------------

    def print_equation_states(self) -> None:
        """ Prints the states of the equations to the screen.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def format_state(state0: tuple) -> str:
            """ Formats the state to be in clear and readable format.

                :param state0: A tuple that represents a state.

                :return: A string representing the state, with its multiplicity.
            """

            # Do not bother with empty states.
            if state0 is None or state0 == tuple():
                return ""

            # Get the state.
            state0_ = state0[0]

            # Get the multiplicity.
            multiplicity0 = state0[1]

            # Format the state.
            state0_ = "<" + ",".join(map(lambda x: f"{x[0]}{x[1]}", state0_)) + ">"

            # Get the multiplicity and attach it if necessary.
            state0_ = str(multiplicity0) + state0_ if multiplicity0 > 1 else state0_

            return state0_

        def get_column_widths(column_names0: dict, equation0: list) -> tuple:
            """ Gets a tuple of integers with the widths of the different
                columns to generate the table.

                :param column_names0: The dictionary with the names of the
                 columns.

                :param equation0: The list of equations for a given state.

                :return: A tuple of integers with the widths of the different
                 columns to generate the table.
            """

            # Get the keys and sort them.
            keys0 = sorted(column_names0.keys(), reverse=True)

            # Get the lengths of the strings.
            lengths0 = list(map(lambda x: len(column_names0[x]), keys0))

            # ---------------------- For the decay states ----------------------

            # Get the strings representing the decay states.
            for key0 in equation0[1].keys():
                # Format all the states of the create states.
                states0 = ", ".join(list(map(lambda x: format_state(x), equation0[2][key0])))

                # Format all the states of the decay states.
                states0_ = ", ".join(list(map(lambda x: format_state(x), equation0[1][key0])))

                # Get the maximum lengths.
                lengths0 = list(map(lambda x, y: max(x, len(y)), lengths0, [key0, states0, states0_]))

            return tuple(lengths0)

        def get_equation(key0: str, equation0: tuple, column_widths0: tuple) -> str:
            """ Returns the string that represents the terms in the equation.

                :param key0: The string representiaton of the particular rate
                 constant for which the rate list is to be obtained.

                :param equation0: The list of equations for a given state.

                :param column_widths0:The widths of each column to properly
                 print the lines and separators.

                :return: The string that represents the terms in the equation.
            """

            # Format all the states of the create states.
            states0 = ", ".join(list(map(lambda x: format_state(x), equation0[2][key0])))

            # Format all the states of the decay states.
            states0_ = ", ".join(list(map(lambda x: format_state(x), equation0[1][key0])))

            # List where the strings will be stored.
            strings0 = [key0, states0, states0_]

            # Format the strings.
            for i0, string0 in enumerate(strings0):
                # If it is the first string.
                if i0 == 0:
                    # Justify the text in the center.
                    strings0[i0] = f"{strings0[i0]:^{column_widths0[i0]}}"

                    continue

                # Justify the text towards the left.
                strings0[i0] = f"{strings0[i0]:<{column_widths0[i0]}}"

            # Join the strings.
            strings0 = get_separator(column_widths) + "\n" + "|" + "|".join(strings0) + "|" + "\n"

            return strings0

        def get_header(column_names0: dict, column_widths0: tuple) -> str:
            """ Returns the formatted header of the table.

                :param column_names0: The dictionary with the names of the
                 columns.

                :param column_widths0:The widths of each column to properly
                 print the lines and separators.

                :return: The string that represents the header of the table.
            """

            # Get the header.
            header0 = sorted(column_names0.keys(), reverse=True)

            # Get the formatted strings.
            header0 = [f"{column_names0[key0]:^{column_widths0[i0]}}" for i0, key0 in enumerate(header0)]

            # Join the strings properly.
            header0 = get_separator(column_widths) + "\n" + "|" + "|".join(header0) + "|"

            return header0

        def get_separator(column_widths0: tuple) -> str:
            """ Gets a string with the row separator.

                :param column_widths0: The widths of each column to properly
                 print the lines and separators.

                :return: A string with the row separator.
            """

            # Get the separator string.
            separator0 = ["-" * length0 for length0 in column_widths0]

            # Join the separators.
            separator0 = "-" + "-".join(separator0) + "-"

            return separator0

        def validate_equations() -> bool:
            """ Validates the equations list is not empty.

                :return: True, if there are equations to display. False,
                 otherwise.
            """

            # Check that the equations list is not empty.
            if len(self.equations) == 0:
                print("There are no equations to print.")

                return False
            
            return True

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        # Validate the equations.
        if not validate_equations():
            return

        # ----------------------------------------------------------------------
        # Define the names of the columns.
        # ----------------------------------------------------------------------

        # The dictionary with the names of the columns.
        column_names = {
            "create": "Create Processes",
            "decay": "Decay Processes",
            "key": "Rate Constant"
        }

        # ----------------------------------------------------------------------
        # Get the column contents.
        # ----------------------------------------------------------------------

        # For every equation in the equations list.
        for i, equation in enumerate(self.equations):

            # # Print the state.
            print("State: ", format_state((equation[0], 1)))

            # Get the column widths.
            column_widths = get_column_widths(column_names, equation)

            # Print the header.
            str_ = get_header(column_names, column_widths) + "\n"

            # For every key.
            for j, key in enumerate(equation[1].keys()):
                # Get the equation terms.
                str_ += get_equation(key, equation, column_widths)

            # Print the separator.
            str_ += get_separator(column_widths) + "\n"

            # Print the table.
            print(str_)

            # Release the memory.
            del str_

    # --------------------------------------------------------------------------
    # Save Methods.
    # --------------------------------------------------------------------------

    @abstractmethod
    def save_equations(self, file_name: str = None, format_type: str = None, order: int = 0, save_path: str = None) -> None:
        """ Generates the equations in the requested format, with the terms
            approximated to the given order.

            :param file_name: The name of the file where the equations are to be
             saved; must be extensionless. Named equations by default.

            :param format_type: A string that represents the format of the
             requested equations. NOT case sensitive, e.g., "A" = "a".

            :param order: An integer that represents the order to which the
             equations must be approximated. Zeroth order, or less, means that
             the equations will be written in an exact way.

            :param save_path: The path where a file with the equations is to be
             created, if at all. None, by default.
        """
        pass

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Private Interface.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Get Methods.
    # --------------------------------------------------------------------------

    @abstractmethod
    def _get_associated_operations(self) -> tuple:
        """ Returns a tuple with the 3-tuples that contain the operations that
            can be applied to the states of the system, the order of the process
            and the string representation of the state.

            :return: A tuple with the 3-tuples that contain the operations that
             can be applied to the states of the system, the order of the
             process and the string representation of the state, in the order:
             (process order, process rate constant, pointer to function).
        """

        # Define the operations.
        process_information = (
            # (process order, process rate constant, process function)
            (None, None, None),
        )

        return process_information

    @abstractmethod
    def _get_constraints(self, states: tuple) -> tuple:
        """ Given a set of numbered states, it gets the constraints of the
            system, i.e., the probability identities, 1 = sum(x in X) P(x),
            P(x) = sum(y in Y) P(x,y), P(y) = sum(x in X) P(x,y), etc; with
            0 <= P(x) <= 1 for x in X.

            :param states: A tuple that containst all of the "left-hand" states
             of the system.

            :return: The constrainsts of the system as equalities.
        """
        pass

    def _get_contracted_state(self, states: list, index: int = -1) -> tuple:
        """ From a list of states, it returns the completely contracted state.
            For this to happen, the list of states must contain as much states
            as there arr possible number of states and the indexes of ALL the
            states must be the same.

            :param states: The list of states that are to be contracted.

            :param index: The index of the list to be contracted. It must be an
             integer number between zero and the length of one of the states, or
             a negative integer number, i.e., an index that indicates the index
             of the array to be contracted.

            :return: A tuple with the contracted state and the original states
             that were contracted; if the return value is an empty tuple, it
             means that the state cannot be contracted.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def get_particles_and_indexes(states0: list, indexes0: list, particles0: list) -> None:
            """ Decomposes the states into individual particles and indexes
                tuples, then appends them to the indexes and particles list.

                :param states0: The list of the states to be contracted.

                :param indexes0: The list in which the indexes must be added.

                :param particles0: The list where the particle of the states are
                 stored.
            """

            # Go through each state.
            for i0, state0 in enumerate(states0):
                # Get the particles and indexes of the state.
                particles0_, indexes0_ = self._get_state_elements(state0)

                # Append them to their respective lists.
                indexes0.append(indexes0_)
                particles0.append(list(particles0_))

        def validate_states(states0: list, index0: int) -> bool:
            """ Checks that the states are valid to perform an index
                contraction operation.

                :param states0: The list of the states to be contracted.

                :param index0: The index of the list to be contracted. It must
                 be a number between zero and the length of one of the states,
                 or a negative number, i.e., an index that indicates the
                 component of the states to be contracted. The length of the
                 states to be contracted must be the same.

                :return: True, if the state is valid for contraction. False,
                 otherwise.
            """

            # All the states are the same length and the number of length of all the states is the same.
            is_valid0 = len(states0) == len(self.states) and len(set(map(lambda x: len(x), states0))) == 1

            # For every state.
            for state0 in states0:
                # Validate the state.
                self._validate_state(state0)

            # Raise an error if the index is not valid.
            if index0 >= len(states0[0]) or not isinstance(index0, int):
                raise ValueError(
                    "The requested entry for contraction must be a negative integer or a positive integer less than"
                    f" {len(states0[0])}. Current value: {index0}, Type: {type(index0)}."
                )

            return is_valid0

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        # While the index is negative.
        while index < 0:
            # Add the length of any of the states.
            index += len(states[0])

        # Validate the states before attempting the contraction.
        if not validate_states(states, index):
            return tuple([]), states

        # ----------------------------------------------------------------------
        # Set the variables where the information will be stored.
        # ----------------------------------------------------------------------

        # Auxiliary variables.
        indexes = []
        particles = []

        # Get the particles and indexes.
        get_particles_and_indexes(states, indexes, particles)

        # ----------------------------------------------------------------------
        # Validate the indexes for contraction.
        # ----------------------------------------------------------------------

        # If there are different indexes do not continue.
        if not len(set(indexes)) == 1:
            return tuple([]), states

        # Get the particles in the column to be contracted.
        particles_ = [states[i][index][0] for i, _ in enumerate(self.states)]

        # Check that the desired row contains all the states.
        if not set(particles_) == set(self.states):
            return tuple([]), states

        # Turn the particles states into a list.
        particles = list(map(list, particles))

        # Remove the given component.
        for j, _ in enumerate(self.states):
            particles[j].pop(index)
        particles = list(map(tuple, particles))

        # Check that the rest of the rest of the components have the same element.
        if not len(set(particles)) == 1:
            return tuple([]), states

        # The contracted index can be the zeroth index.
        indexes = list(indexes[0])

        # Without the given index.
        indexes.pop(index)

        # Remember to use the probability identity.
        if len(indexes) == 0:
            return (1,), states

        # Get the contracted state.
        contracted = tuple((state_, indexes[i]) for i, state_ in enumerate(particles[0]))

        return contracted, states

    def _get_decay_states(self, state: tuple, operations: dict) -> tuple:
        """ Given a state and a set of operations, it returns the states that
            are a generated when ALL the possible operations are performed on
            a given state.

            :param state: A tuple that represents a state on which to operate.

            :param operations: A dictionary with all the possible operations of
             the system.

            :return: A 2-tuple that contains the original state and a dictionary
             with all the possible states generated by the specific collection
             of operations.
        """

        # Validate the state.
        self._validate_state(state)

        # ----------------------------------------------------------------------
        # Initialize the variables.
        # ----------------------------------------------------------------------

        # Initialize a dictionary.
        decay_states = {}

        # Get the dictionary keys.
        keys = operations.keys()

        # ----------------------------------------------------------------------
        # Get the decay states.
        # ----------------------------------------------------------------------

        # For each process.
        for key in keys:
            # Get all the decay states.
            decay_states[key] = operations[key](state)

        return state, decay_states

    def _get_is_substate(self, state1: tuple, state2: tuple) -> bool:
        """ Determines if state 1 is a substate, or proper substate, of state 2;
            they must be in the same format; order matters in this case.

            :param state1: The state that is to be found within state2.

            :param state2: The state where state1 is going to be searched.

            :return: True if state1 is a substate, or proper substate, of
             state2. False, otherwise.
        """

        # Validate the first state.
        self._validate_state(state1)

        # Validate the second state.
        self._validate_state(state2)

        # State 1 cannot be a substate of state 2.
        if len(state1) > len(state2):
            return False

        # If the states are equal no need to continue.
        if state1 == state2:
            return True

        # Sort substate 1 in index order.
        substate1 = sorted(state1, key=lambda x: (x[1], x[0]))

        # Get all the possible combinations of the components in state 2.
        substate2 = list(map(tuple, itertools.combinations(state2, len(substate1))))

        # For every combination of state 2.
        for substate in substate2:
            # If an equality is found.
            if substate1 == sorted(substate, key=lambda x: (x[1], x[0])):
                return True

        return False

    def _get_multiplicity(self, state_dictionary: dict) -> dict:
        """ Given the decay/created states dictionary, it returns a dictionary
            with the UNIQUE states for each process and their multiplicity.

            :param state_dictionary: The dictionary of processes associated with
             a state.

            :return: A dictionary with the UNIQUE states for each process and
             their multiplicity.
        """

        # Auxiliary variables.
        keys = state_dictionary.keys()

        # The dictionary that contains the UNIQUE processes and multiplicities.
        multiplicity_dictionary = {}

        # For each process.
        for key in keys:
            # Get the unique states.
            unique_states = set(state_dictionary[key])

            # Validate the states.
            for state in state_dictionary[key]:
                self._validate_state(state)

            # Get the new dictionary entry.
            multiplicity_dictionary[key] = [(state, state_dictionary[key].count(state)) for state in unique_states]

        # Return the multiplicity dictionary.
        return multiplicity_dictionary

    @abstractmethod
    def _get_numbering(self, state):
        """ Returns a tuple with the possible numbering a state of length N can
            have. The format of a SINGLE state must be given in the format:

            ( (particle_at_site1, numbering_scheme1),
                            .
                            .
                            .
              (particle_at_siteN, numbering_schemeN),
            )

            :param state: The state to be numbered.

            :return: The list of possible numbered states in the given format.
        """
        pass

    @abstractmethod
    def _get_process_functions(self):
        """ Returns all the pointers to the functions that operate on the
            different states to, potentially, modify them.

            :return process_functions: A dictionary, whose keys are the string
            represenation of the rates, with the functions that will potentially
            modify a state.
        """

        # Get the associated operations.
        operations = self._get_associated_operations()

        # Get the rate strings tuple.
        process_functions = tuple((operation[1], operation[2],) for operation in operations)

        # Convert it into a dictionary.
        process_functions = dict(process_functions)

        return process_functions

    @abstractmethod
    def _get_process_orders(self):
        """ Returns all the dictionary of integers that represent the minimum
            number of sites required for the given processes to take place.

            :return rates_order: A dictionary, whose keys are the string
            represenation of the rates, with the integers that represent the
            minimum number of sites required for the given processes to take
            place.
        """

        # Get the associated operations.
        operations = self._get_associated_operations()

        # Get the rate strings tuple.
        rates_orders = tuple((operation[1], operation[0], ) for operation in operations)

        # Convert it into a dictionary.
        rates_orders = dict(rates_orders)

        return rates_orders

    @abstractmethod
    def _get_process_rates(self):
        """ Returns all the string representation of the rates associated with
            the class.

            :return rates_strings: A tuple with the string rates associated with
            the class.
        """

        # Get the associated operations.
        operations = self._get_associated_operations()

        # Get the rate strings tuple.
        rates_strings = tuple(operation[1] for operation in operations)

        return rates_strings

    def _get_products_create(self, state, states_decay):
        """ Given a state and the list of states, whose decay states due to the
            different process must be included, it returns of the states that
            create the given state.

            :param state: The state that must be created from the other states.

            :param states_decay: A list of states and their decay products. It
            must be in the format.
                (state, 2-tuple with state and dictionary with decay states)

            :return: A 2-tuple with the state and a dictionary with the UNIQUE
            lowest order states that will create the state through a given
            process.
        """

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Get the keys to the dictionary.
        keys = states_decay[0][1].keys()

        # Start an empty dictionary of lists, from the keys.
        create_dictionary = {key: [] for key in keys}

        # For all the decay states.
        for state_decay in states_decay:
            # Only non-substates are possible..
            if self._get_is_substate(state, state_decay[0]):
                continue

            # For all the processes.
            for key in keys:
                # For all the states formed by a particular process.
                for state_0 in state_decay[1][key]:
                    # Only possible if the state appears in the original.
                    if not self._get_is_substate(state, state_0):
                        continue

                    # Add the state to the dictionary.
                    create_dictionary[key].append(state_decay[0])

        # Reduce the entries of the dictionary.
        create_dictionary = self._reduce_to_unique_states(create_dictionary, state)

        # Reduce the entries of the dictionary.
        create_dictionary = self._get_multiplicity(create_dictionary)

        return create_dictionary

    def _get_products_decay(self, state, states_decay):
        """ Given a state and the list of states, whose decay states due to the
            different process must be included, it returns of the states that
            make the first state decay.

            :param state: The state whose decay process are to be obtained.

            :param states_decay: A list of states and their decay products. It
            must be in the format.
                (state, 2-tuple with state and dictionary with decay states)

            :return: A 2-tuple with the state and a dictionary with the UNIQUE
            lowest order states that will make it decay through a given process.
        """

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Get the keys to the dictionary.
        keys = states_decay[0][1].keys()

        # Start an empty dictionary of lists, from the keys.
        decay_dictionary = {key: [] for key in keys}

        # For all the decay states.
        for state_decay in states_decay:
            # Only substates are possible..
            if not self._get_is_substate(state, state_decay[0]):
                continue

            # For all the processes.
            for key in keys:
                # For all the states formed by a particular process.
                for state_0 in state_decay[1][key]:
                    # Only if the state does not appear in the original.
                    if self._get_is_substate(state, state_0):
                        continue

                    # Add the state to the dictionary.
                    decay_dictionary[key].append(state_decay[0])

        # Reduce the entries of the dictionary.
        decay_dictionary = self._reduce_to_unique_states(decay_dictionary, state)

        # Reduce the entries of the dictionary.
        decay_dictionary = self._get_multiplicity(decay_dictionary)

        return decay_dictionary

    def _get_states(self, order: int = 1) -> tuple:
        """ Given the order, it returns a list of ALL the possible combinations
            of the system variables, i.e., all the possible combinations of the
            variables in N slots, where N=order; NON-NUMBERED.

            :param order: The order of the requested states.

            :return: A list of all the possible states of the given order.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def validate_order():
            """ Validates that the order to which the equations are requested
                is valid.
            """

            # Check the requested order is greater than zero.
            if order <= 0:
                raise ValueError(
                    "The order parameter must be greater than zero."
                )

            # Check that the order parameter is not more than the number of sites.
            if order > self.sites_number:
                raise ValueError(
                    f"The order parameter must less than or equal to {self.sites_number}."
                )

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Validate that the order is valid.
        validate_order()

        # Get an iterator to get the states.
        all_states = product(*[cp.deepcopy(self.states) for _ in range(order)])

        # Make the states into a list.
        all_states = list(all_states)

        return all_states

    def _get_states_left(self, order):
        """ Returns a list of the numbered lowest order states for the
            equations to be written in; i.e., the states that have the
            differential operator d/dt.

            :param order: The order to which the equation is to be written.

            :return states: A list of the non-numbered lowest order states for
            the equations to be written in.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def validate_order():
            """ Validates that the order is an integer greater than or equal to
                zero.
            """

            # If the order is not valid.
            if order < 0 or not isinstance(order, (int,)):
                raise ValueError(
                    f"The order must be an integer number greater than or equal"
                    f" to zero. Current order = {order}, Type: {type(order)}"
                )

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Validate the order.
        validate_order()

        # Auxiliary variables.
        states = []
        states_0 = []

        # Fix the order, if needed.
        order0 = 1 if order == 0 else order
        order0 = order0 if order < self.sites_number else self.sites_number

        # Get the states up to the given order.
        if order0 < self.sites_number:
            # Up until the requested order.
            for i in range(1, order0 + 1):
                # Get the un-numbered states.
                states_0.extend(self._get_states(i))
        else:
            # Exact Equations.
            states_0.extend(self._get_states(self.sites_number))

        # For every un-numbered state.
        for state_0 in states_0:
            # Get the numbered states.
            states.extend(self._get_numbering(state_0))

        return states

    def _get_states_right(self, order):
        """ Returns a list of the numbered states that will potentially appear
            in the derivative term.

            :param order: The order to which the equation is to be written.

            :return states: A list of the non-numbered lowest order states for
            the equations to be written in.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def validate_order():
            """ Validates that the order is an integer greater than or equal to
                zero.
            """

            # If the order is not valid.
            if order < 0 or not isinstance(order, (int,)):
                raise ValueError(
                    f"The order must be a number greater than or equal to zero."
                    f" Current order = {order}, Type = {type(order)}"
                )

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Validate the order.
        validate_order()

        # Auxiliary variables.
        orders = []
        states = []
        states_0 = []

        # Fix the order, if needed.
        order0 = 1 if order == 0 else order
        order0 = order0 if order < self.sites_number else self.sites_number

        # Get the states up to the given order.
        if order0 < self.sites_number:
            # Get the orders of the process.
            processes_orders = list(set(order_0[0] for order_0 in self._get_associated_operations()))

            # For each length of states.
            for order_0 in range(1, order0 + 1):
                # Get the maximum order.
                orders.extend([order_0 + order_1 - 1 for order_1 in processes_orders])

            # Order vectors cannot be longer than the total number of sites.
            maximum_order = max(orders)
            maximum_order = min(maximum_order, self.sites_number)

            # For all the orders, up to the maximum order.
            for i in range(1, maximum_order + 1):
                # Get the states.
                states_0.extend(self._get_states(i))

        else:
            # Exact Equations.
            states_0.extend(self._get_states(self.sites_number))

        # For every non-numbered state
        for state0 in states_0:
            # Number the states associated with the non-numbered state.
            states.extend(self._get_numbering(state0))

        return states

    def _get_state_elements(self, state):
        """ Gets two tuples one with the elements of the given state; one with
            all the particles and the other one with the indexes.

            :param state: The representation of the state, that must be in
            the standard form.

            :return: A 2-tuple that contains the particles in the state and its
            numerical indexes in the order:

            (particles, numerical indexes)
        """

        # Validate the state.
        self._validate_state(state)

        # Get the numpy representation of the particles of the state.
        particle_list = tuple(entry[0] for entry in state)

        # Get the numpy representation of the indexes of the state.
        index_list = tuple(entry[1] for entry in state)

        return particle_list, index_list

    # --------------------------------------------------------------------------
    # Other Methods.
    # --------------------------------------------------------------------------

    @abstractmethod
    def _reduce_to_unique_states(self, state_list, target_state):
        """ Given a list of states, it attempts to contract

            :param state_list: The list of states to be reduced.

            :param target_state: The state that is being targeted to appear in
            the reduced list.

            :return: A list of the reduced states in the format (state,
            multiplicity).
        """
        pass

    # --------------------------------------------------------------------------
    # Validation Methods.
    # --------------------------------------------------------------------------

    def _validate_sites_number(self) -> None:
        """ Validates that the constrainst of the system have the proper format.
        """

        # Verify the number of sites is a positive integer.
        if self.sites_number < 1:
            raise ValueError(
                "The number of sites must an integer greater than zero."
            )

    def _validate_state(self, state: tuple) -> None:
        """ Validates that the state is a collection of 2-tuples and that
            the tuples have the form ("particle", "site"); where "particle" is
            in the set self.states and "site" is an integer in the range
            [1, self.sites].

            :param state: The state to validate.
        """

        # Make a copy of the state.
        state = cp.deepcopy(state)

        # ----------------------------------------------------------------------
        # Check that the state is made of tuples of length 2.
        # ----------------------------------------------------------------------

        # Check that each entry in the state is a tuple.
        if not all(map(lambda x: isinstance(x, (tuple,)), state)):
            tmp_types = [str(type(x)) for x in state]
            raise TypeError(
                f"A state must be a collection of tuples, at least one elements"
                f" is not a tuple. State: {state}, Types of state: {tmp_types}."
            )

        # Check that each entry in the state is 2 sites long.
        if not all(map(lambda x: len(x) == 2, state)):
            tmp_lengths = [str(len(x)) for x in state]
            raise TypeError(
                f"All states must be tuples of length 2. Tuple lengths:"
                f" {tmp_lengths}"
            )

        # Check that the maximum length of the state is the number of sites.
        if not 1 <= len(state) <= self.sites_number:
            raise ValueError(
                f"The current length of the state list is not valid, it must be"
                f"in the range [1, {self.sites_number}].  Current legth:"
                f" {len(state)}."
            )

        # ----------------------------------------------------------------------
        # Set the auxiliary variables.
        # ----------------------------------------------------------------------

        # Get the particles.
        state = np.array(state)

        # Get the particles and the numbering of the states.
        particles = list(state[:, 0])
        numbering = tuple([int(x[1]) for x in state])

        # ----------------------------------------------------------------------
        # Check the state.
        # ----------------------------------------------------------------------

        # Check that the states are made of valid particles.
        if not all(map(lambda x: x in self.states, particles)):
            raise ValueError(
                f"The states are not valid, they must be in the list "
                f"{self.states}. Current particles in the state: {particles}."
            )

        # ----------------------------------------------------------------------
        # Check the numbering.
        # ----------------------------------------------------------------------

        # Check that the length of the state is greater than zero and less than or equal to the maximum number of sites.
        if not 1 <= len(numbering) <= self.sites_number:
            raise ValueError(
                f"The current length of the numbering list is not valid, it "
                f"must be in the range [1, {self.sites_number}].  Current "
                f"length: {len(numbering)}."
            )

        # Check that the numbering of the state is greater than zero and less than or equal to the maximum number of
        # sites.
        if not all(map(lambda x: 0 < x <= self.sites_number, numbering)):
            raise ValueError(
                f"The numbering of the states must be in the range [1,"
                f" {self.sites_number}]. There is an index that is not in this "
                f"range: {numbering}."
            )

        # Check that the numbering for each site is unique.
        if not len(set(numbering)) == len(numbering):
            raise ValueError(
                f"Indexes in the numbering list MUST be unique. There is a "
                f"non-unique index: {numbering}."
            )

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Constructor, Dunder Methods and Dunder Variables.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # -------------------------------------------------------------------
    # Constructor.
    # --------------------------------------------------------------------------

    def __init__(self, sites_number: int, states: Union[list, tuple]):
        """ Creates an EquationGenerator object and initializes its properties.

            :param sites_number: The integer number of sites that the system
             has.

            :param states: A list of unique strings, that represent the names of
             the statistically independent variables that each site of the
             system can take.
        """

        # ----------------------------------------------------------------------
        # Define the default model parameters.
        # ----------------------------------------------------------------------

        # Define the number of sites the system has.
        self.sites_number = sites_number

        # Define the possible unique states each site of the system can take.
        self.states = states

        # Define the list where the constraint equations are saved.
        self.constraints = []

        # Define the list where the equations are saved.
        self.equations = []

    # --------------------------------------------------------------------------
    # Dunder Methods.
    # --------------------------------------------------------------------------
