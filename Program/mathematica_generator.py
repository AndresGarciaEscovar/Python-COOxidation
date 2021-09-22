"""
"""



# Imports.
import copy as cp
import itertools
import numpy as np

from abc import ABC, abstractmethod
from collections.abc import Iterable
from itertools import product


class EquationGenerator(ABC):
    """ Generates the differential equations in different formats. Currently
        only Mathematica and LaTeX are supported.

        :param self.equations: The list where the equations will be saved.

        :param self.sites: The maximum number of sites

        :param  self.states: The UNIQUE states in which each side can be in.
    """

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Getters, Setters and Deleters.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    @property
    def number_of_sites(self):
        """ Gets the number of sites in the lattice.
        """
        return cp.deepcopy(self.__number_of_sites)

    @number_of_sites.setter
    def number_of_sites(self, number_of_sites):
        """ Sets the number of sites in the lattice.
        """

        # Verify that the parameter is a positive integer number.
        if not isinstance(number_of_sites, (int,)) or number_of_sites < 1:
            raise ValueError("The number of sites must be an integer greater than zero.")

        self.__number_of_sites = number_of_sites

    @number_of_sites.deleter
    def number_of_sites(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def states(self):
        """ Gets the possible states a system can have.
        """
        return cp.deepcopy(self.__states)

    @states.setter
    def states(self, states):
        """ Sets the possible states a system can have. Must be an iterable
            object of UNIQUE objects that allow a string representation.
        """

        # Verify that the states variable can be iterated through.
        if not isinstance(states, Iterable):
            raise ValueError(f"The states variable must be an iterable object, current type: {type(states)}")

        # Verify the states allow a string representation.
        tmp_states = tuple(map(str, states))

        # Check that the elements are unique.
        if not len(tmp_states) == len(set(states)):
            raise ValueError(f"The states a system can take must be unique: {states}")

        self.__states = tmp_states

    @states.deleter
    def states(self):
        pass

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Public Interface.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Generate methods.
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # Get methods.
    # --------------------------------------------------------------------------

    def get_nth_order_equations(self, processes, order=0, print_equations=False):
        """ Gets the nth order equations for the system in the format
                (state, decay_states, create_states),
            where "decay_states" and "create_states" are the dictionaries that
            contain the states to which the "state" decays, or the states from
            where the state is created, due to the different process in the
            system.

            :param order: The lowest order to which the equations must be
            given.

            :param processes: The information of the processes. This must be
            a tuple of tuples, such that the inner tuples are in the format
            (process order, process rate constant, pointer to function).

            :param print_equations: If the table of equations must be printed.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Get the lowest order states.
        states_left_hand = self._get_states_left(order)

        # Get the other involved states.
        states_right_hand = self._get_states_right(order)

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Private Interface.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Generate methods.
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # Get methods.
    # --------------------------------------------------------------------------

    @abstractmethod
    def _get_associated_operations(self):
        """ Returns a tuple with the 3-tuples that contain the operations that
            can be applied to the states of the system, the order of the process
            and the string representation of the state.

            :return process_information: A tuple with the 3-tuples that contain
            the operations that can be applied to the states of the system, the
            order of the process and the string representation of the state, in
            the order:

            (process order, process rate constant, pointer to function)
        """

        # Define the operations.
        process_information = (
            # (process order, process rate constant, process function)
            (None, None, None),
        )

        return process_information

    def _get_contracted_state(self, states, entry=-1):
        """ From a list of states, it returns the completely contracted state.
            For this to happen, the list of states must contain as much states
            as there arr possible number of states and the indexes of ALL the
            states must be the same.

            :param states: The list of states that are to be contracted.

            :param entry: The entry of the list to be contracted. It must be an
            integer number between zero and the length of one of the states, or
            a negative integer number, i.e., an index that indicates the index
            of the array to contracted.

            :return: A tuple with the contracted state and the original states
            that were contracted; if the return value is an empty tuple, it
            means that the state cannot be contracted.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def get_states_and_indexes():
            """ Decomposes the states into individual particles and indexes
                tuples, then appends them to the indexes and particles list.
            """

            # Go through each state.
            for i, state in enumerate(states):
                # Get the particles and indexes of the state.
                tmp_particles0, tmp_indexes0 = self._get_state_elements(state)

                # Append them to their respective lists.
                states_indexes.append(tmp_indexes0)
                states_particles.append(list(tmp_particles0))

        def validate_states(entry0):
            """ Checks that the states are valid to perform an index
                contraction operation.

                :param entry0: The entry0 of the list to be contracted. It must
                be a number between zero and the length of one of the states, or
                a negative number, i.e., an index that indicates the index of
                the array to contracted. The length of the states to be
                contracted must be the same.

                :return True, if the state is valid for contraction. False,
                otherwise.
            """

            # Verify it is a list of valid states.
            for state0 in states:
                self._validate_state(state0)

            # Verify that the list contains as much states as possible site states.
            if not len(states) == len(self.states):
                raise ValueError("When a state is to be contracted, there must be as many states as there are"
                                 f" particles. Number of requested = {len(states)}, possible site states ="
                                 f" {len(self.states)}")

            # Verify that the index is within the limits.
            while entry0 < 0 and isinstance(entry, int):
                entry0 += len(states[0])

            # Raise an error if the entry is not valid.
            if entry0 >= len(states[0]) or not isinstance(entry, int):
                raise ValueError("The requested entry for contraction must be a negative integer or a positive"
                                 f" integer less than {len(states[0])}. Current value: {entry}, Type: {type(entry)}.")

            # ----------------------------------------------------------------------
            # Return if the index is valid.
            # ----------------------------------------------------------------------

            # All the states must have the same length to continue the process.
            if not len(set(map(lambda x: len(x), states))) == 1:
                return False

            return True

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Validate the states
        if not validate_states(entry):
            return tuple([]), states

        # Auxiliary variables.
        states_indexes = []
        states_particles = []

        # Get the states and indexes.
        get_states_and_indexes()

        # ----------------------------------------------------------------------
        # Validate the indexes for contraction.
        # ----------------------------------------------------------------------

        # If there are different indexes do not continue.
        if not len(set(states_indexes)) == 1:
            return tuple([]), states

        # Get the row to be contracted.
        tmp_particles = [states[i][entry][0] for i, _ in enumerate(self.states)]

        # Check that the desired row contains all the states.
        if not set(tmp_particles) == set(self.states):
            return tuple([]), states

        # Turn the particles states into a list.
        states_particles = list(map(list, states_particles))

        # Remove the given entry.
        for j, _ in enumerate(self.states):
            states_particles[j].pop(entry)
        states_particles = list(map(tuple, states_particles))

        # Check that the rest of the rest of the entries have the same element.
        if not len(set(states_particles)) == 1:
            return tuple([]), states

        # Remove the contracted index.
        states_indexes = list(states_indexes[0])
        states_indexes.pop(entry)

        # Remember to use the probability identity.
        if len(states_indexes) == 0:
            return (1,), states

        # Get the contracted state.
        contracted_state = tuple((state_0, states_indexes[i]) for i, state_0 in enumerate(states_particles[0]))

        return contracted_state, states

    def _get_decay_states(self, state, operations):
        """ Given a state and a set of operations, it returns the states that
            are a generated when ALL the possible operations are performed on
            a given state.

            :param state: The state on which to operate.

            :param operations: A dictionary with all the possible operations of
            the system.

            :return: A 2-tuple that contains the original state and a dictionary
            with all the possible states generated by the specific collection of
            operations.
        """

        # Validate the state.
        self._validate_state(state)

        # Auxiliary variables.
        keys = operations.keys()

        # Initialize a dictionary.
        decay_states = {}

        # For each process.
        for key in keys:
            # Get all the decay states.
            decay_states[key] = operations[key](state)

        return state, decay_states

    def _get_is_substate(self, state1, state2):
        """ Determines if state 1 is a substate, or proper substate, of state 2;
            they must be in the same format; order matters in this case.

            :param state1: The state that is to be found within state2.

            :param state2: The state where state1 is going to be searched.

            :return: True if state1 is a substate, or proper substate, of
            state2. False, otherwise.
        """

        # Validate the states.
        self._validate_state(state1)
        self._validate_state(state2)

        # State 1 cannot be a substate of state 2.
        if len(state1) > len(state2):
            return False

        # If the states are equal no need to continue.
        if state1 == state2:
            return True

        substate1 = sorted(state1, key=lambda x: (x[1], x[0]))
        substate2 = list(map(tuple, itertools.combinations(state2, len(substate1))))

        # Otherwise, ALL the entries in state1 must be in state2.
        for substate in substate2:
            # If an equality is found.
            if substate1 == sorted(substate, key=lambda x: (x[1], x[0])):
                return True

        return False

    def _get_multiplicity(self, state_dictionary):
        """ Given the decay/created states dictionary, it returns a dictionary
            with the UNIQUE states for each process and their multiplicity.

            :param state_dictionary: The dictionary of processes associated with
            a state.

            :return multiplicity_dictionary: A dictionary with the UNIQUE states
            for each process and their multiplicity.
        """

        # Auxiliary variables.
        keys = state_dictionary.keys()

        # Define the dictionary that contains the UNIQUE processes and their
        # multiplicity.
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

    def _get_states(self, order=1):
        """ Given the order, it returns a list of ALL the possible combinations
            of the system variables, i.e., all the possible combinations of the
            variables in N slots, where N=order; NON-NUMBERED.

            :param order: The order of the requested states.

            :return all_states: A list of all the possible states of the given
            order.
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
                raise ValueError("The order parameter must be greater than zero.")

            # Check that the order parameter is not more than the number of sites.
            if order > self.number_of_sites:
                raise ValueError(f"The order parameter must less than or equal to {self.number_of_sites}.")

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
                raise ValueError(f"The order must be an integer number greater than or equal to zero."
                                 f" Current order = {order}, Type: {type(order)}")

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
        order0 = order0 if order < self.number_of_sites else self.number_of_sites

        # Up until the requested order.
        for i in range(1, order0 + 1):
            # Get the un-numbered states.
            states_0.extend(self._get_states(i))

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
                raise ValueError(f"The order must be a number greater than or equal to zero."
                                 f" Current order = {order}, Type = {type(order)}")

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Validate the order.
        validate_order()

        # Auxiliary variables.
        orders = []
        states = []

        # Fix the order, if needed.
        order0 = 1 if order == 0 else order
        order0 = order0 if order < self.number_of_sites else self.number_of_sites

        # Get the orders of the process.
        processes_orders = list(set(order_0[0] for order_0 in self._get_associated_operations()))

        # For each length of states.
        for order_0 in range(1, order0 + 1):
            orders.extend([order_0 + order_1 - 1 for order_1 in processes_orders])

        # Order vectors cannot be longer than the total number of sites.
        maximum_order = max(orders)
        maximum_order = min(maximum_order, self.number_of_sites)

        # Get the states of the maximum order.
        states_0 = self._get_states(maximum_order)

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
    # Validation methods.
    # --------------------------------------------------------------------------

    def _validate_state(self, state):
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
            raise TypeError(f"A state must be a collection of tuples, at least one elements is not "
                            f"a tuple. State: {state}, Types of state: {tmp_types}")

        # Check that each entry in the state is 2 sites long.
        if not all(map(lambda x: len(x) == 2, state)):
            tmp_lengths = [str(len(x)) for x in state]
            raise TypeError(f"All states must be tuples of length 2. Tuple lengths: {tmp_lengths}")

        # Check that the maximum length of the state is the number of sites.
        if not 1 <= len(state) <= self.number_of_sites:
            raise ValueError(f"The current length of the state list is not valid, it must be"
                             f"in the range [1, {self.number_of_sites}].  Current legth: {len(state)}.")

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
            raise ValueError(f"The states are not valid, they must be in the list {self.states}."
                             f" Current particles in the state: {particles}.")

        # ----------------------------------------------------------------------
        # Check the numbering.
        # ----------------------------------------------------------------------

        # Check that the length of the state is greater than zero and less than or equal to the maximum number of sites.
        if not 1 <= len(numbering) <= self.number_of_sites:
            raise ValueError(f"The current length of the numbering list is not valid, it must be"
                             f"in the range [1, {self.number_of_sites}].  Current length: {len(numbering)}.")

        # Check that the numbering of the state is greater than zero and less than or equal to the maximum number of
        # sites.
        if not all(map(lambda x: 0 < x <= self.number_of_sites, numbering)):
            raise ValueError(f"The numbering of the states must be in the range [1, {self.number_of_sites}]."
                             f" There is an index that is not in this range: {numbering}")

        # Check that the numbering for each site is unique.
        if not len(set(numbering)) == len(numbering):
            raise ValueError(f"Indexes in the numbering list MUST be unique."
                             f" There is a non-unique index: {numbering}")

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Constructor, Dunder Methods and Dunder Variables.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # -------------------------------------------------------------------
    # Constructor.
    # --------------------------------------------------------------------------

    def __init__(self, number_of_sites, states):
        """ Creates an EquationGenerator object and initializes its properties.

            :param number_of_sites: The number of sites that the system has.

            :param states: A list of unique strings, that represent the names of
            the statistically independent variables that each site of the system
            can take.
        """

        # ----------------------------------------------------------------------
        # Define the default model parameters.
        # ----------------------------------------------------------------------

        # Define the maximum number of sites.
        self.number_of_sites = number_of_sites

        # Define the possible unique states each site of the system can take.
        self.states = states

        # Array where the equations are saved.
        self.equations = []

    # --------------------------------------------------------------------------
    # Dunder Methods.
    # --------------------------------------------------------------------------


class EquationFormatter:
    """ A static class that contains equation formatting functions that are
        generated by the equation generator.
    """

    # --------------------------------------------------------------------------
    # CONSTANTS.
    # --------------------------------------------------------------------------

    CONSTANTS = {"latex": "latex", "mathematica": "mathematica"}

    # --------------------------------------------------------------------------
    # Get methods.
    # --------------------------------------------------------------------------

    @staticmethod
    def get_formatter(format_string="latex"):
        """ Gets the pointer to the function of the requested format.

            :param format_string: The string with the name of the requested
            formatter. NOT case-sensitive.

            :return formatter: The pointer to the function of the formatter.
        """

        # Verify that the function exists.
        if format_string.lower() not in EquationFormatter.CONSTANTS.keys():
            raise ValueError("The requested constants must be in the dictionary: ")

        # Return the function that sets the equation in LaTeX format.
        if format_string.lower() == EquationFormatter.CONSTANTS["latex"]:
            return EquationFormatter.get_latex_format

        # Return the function that sets the equation in Mathematica format.
        if format_string.lower() == EquationFormatter.CONSTANTS["mathematica"]:
            return EquationFormatter.get_mathematica_format

    @staticmethod
    def get_latex_format(state, derivative=False, order=0, rate_constant=False):
        """ Gets the string that represents a state in LaTeX format.

            :param state: The state for which to get LaTeX format string. It
            must be an N-tuple of 2-tuples, where the first index represents
            the particle and the second index represents the site where the
            particle is sitting.

            :param derivative: True, if the state must be returned in derivative
            format. False, otherwise.

            :param order: The order to which the state must be expanded. Order
            zero means the state must not be modified. Other orders means the
            state must be mean-field expanded to the given order.

            :param rate_constant: True, if it is the rate constant that must be
            formatted. False, otherwise.

            :return latex_state: The string that represents the state in a
            LaTeX format.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def get_enclosed_state(state0):
            """ Given an N-tuple of 2-tuples, returns the given state in a LaTeX
                format.

                :param state0: An N-tuple of 2-tuples, where the first index
                of the 2-tuples represents the particle and the second index
                represents the site where the particle is sitting.

                :return: A string representation of the state.
            """

            # Enclose the state properly.
            state_string0 = "\\left<" + ",".join([single_term(term) for term in state0]) + "\\right>"

            return state_string0

        def get_nth_order():
            """ Gets the list a states, such that the states are n sites long.
                If the length of the initial state is less than the order of the
                requested approximation, it returns the state itself.

                :return state_list0: A list of states, each n-sites long or
                shorter.
            """

            # Auxiliary variables.
            state_list0 = []

            # Dont approximate the states that are shorter.
            if len(state) <= order:
                return [state]

            # Split the state in the required number of states.
            for i, _ in enumerate(state):
                # Do this process until all the possible states are obtained.
                if i + order > len(state):
                    break

                # Append the states.
                state_list0.append(cp.deepcopy(state[i: i + order]))

            return state_list0

        def get_nth_order_string():
            """ Gets the string for the states expanded to the nth order.

                :return: The formatted string to the nth order.
            """

            # Auxiliary variables.
            numerator_states0 = []
            denominator_states0 = []

            # Run along all the states.
            for i, state0 in enumerate(state_list):

                # Get the numerator states.
                numerator_states0.append(state0)

                # Only need to get this terms for contiguous states
                if i > len(state_list) - 2:
                    continue

                # Get the intersection of the ith state and the (i+1)th state.
                denominator_states0.append(tuple(sorted(tuple(set(state_list[i]) & set(state_list[i + 1])), key=lambda x: x[1])))

            # Get the LaTeX representation for the states.
            numerator_states0 = "".join(list(map(get_enclosed_state, numerator_states0)))
            denominator_states0 = "".join(list(map(get_enclosed_state, denominator_states0)))

            # Create the string.
            state_string0 = "\\frac{" + numerator_states0 + "}{" + denominator_states0 + "}" if order > 1 else numerator_states0

            return state_string0

        def get_rate_constant():
            """ Formats the rate constant such that it matches the LaTex format.

                :return rate_string0: The formatted rate constant in LaTeX form.
            """

            # The rate constant list must be a string.
            rate_constant_list0 = state.split(".")

            # Sub-index the entries.
            closing_string0 = "}" * (len(rate_constant_list0)-1)
            rate_string0 = "_{".join(rate_constant_list0) + closing_string0

            return rate_string0

        def single_term(term: tuple) -> str:
            """ Gets the string representation of a single term.

                :param term: A 2-tuple, where the first index represents the
                particle and the second index represents the site where the
                particle is sitting.

                :return: A string representation of the term.
            """

            # Format the term properly.
            term_string = f"{term[0]}" + "_{" f"{term[1]}" + "}"

            return term_string

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Auxiliary variables.
        state_string = ""

        # If it is the rate constant the one that is to be formatted.
        if rate_constant:
            return get_rate_constant()

        if order == 0 or len(state) <= order:
            # Format the string.
            state_string += get_enclosed_state(state)

        else:
            # Format the string.
            state_list = get_nth_order()
            state_string += get_nth_order_string()

        # Format the state in derivative form.
        if derivative:
            state_string = "\\frac{d" + state_string + "}{dt}"

        return state_string

    @staticmethod
    def get_mathematica_format(state, derivative=False, order=0, rate_constant=False):
        """ Gets the string that represents a state in Mathematica format.

            :param state: The state for which to get Mathematica format string.
            It must be an N-tuple of 2-tuples, where the first index represents
            the particle and the second index represents the site where the
            particle is sitting.

            :param derivative: True, if the state must be returned in derivative
            format. False, otherwise.

            :param order: The order to which the state must be expanded. Order
            zero means the state must not be modified. Other orders means the
            state must be mean-field expanded to the given order.

            :param rate_constant: True, if it is the rate constant that must be
            formatted. False, otherwise.

            :return latex_state: The string that represents the state in a
            Mathematica format.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def get_enclosed_state(state0):
            """ Given an N-tuple of 2-tuples, returns the given state in a
                Mathematica format.

                :param state0: An N-tuple of 2-tuples, where the first index
                of the 2-tuples represents the particle and the second index
                represents the site where the particle is sitting.

                :return: A string representation of the state.
            """

            # Enclose the state properly.
            state_string0 = "P" + "".join([single_term(term) for term in state0]) + "[t]"

            return state_string0

        def get_nth_order():
            """ Gets the list a states, such that the states are n sites long.
                If the length of the initial state is less than the order of the
                requested approximation, it returns the state itself.

                :return state_list0: A list of states, each n-sites long or
                shorter.
            """

            # Auxiliary variables.
            state_list0 = []

            # Dont approximate the states that are shorter.
            if len(state) <= order:
                return [state]

            # Split the state in the required number of states.
            for i, _ in enumerate(state):
                # Do this process until all the possible states are obtained.
                if i + order > len(state):
                    break

                # Append the states.
                state_list0.append(cp.deepcopy(state[i: i + order]))

            return state_list0

        def get_nth_order_string():
            """ Gets the string for the states expanded to the nth order.

                :return: The formatted string to the nth order.
            """

            # Auxiliary variables.
            numerator_states0 = []
            denominator_states0 = []

            # Run along all the states.
            for i, state0 in enumerate(state_list):

                # Get the numerator states.
                numerator_states0.append(state0)

                # Only need to get this terms for contiguous states
                if i > len(state_list) - 2:
                    continue

                # Get the intersection of the ith state and the (i+1)th state.
                denominator_states0.append(tuple(sorted(tuple(set(state_list[i]) & set(state_list[i + 1])), key=lambda x: x[1])))

            # Get the LaTeX representation for the states.
            numerator_states0 = "*".join(list(map(get_enclosed_state, numerator_states0)))
            denominator_states0 = "*".join(list(map(get_enclosed_state, denominator_states0)))

            # Create the string.
            state_string0 = f"{numerator_states0}/({denominator_states0})" if order > 1 else numerator_states0

            return state_string0

        def get_rate_constant():
            """ Formats the rate constant such that it matches the LaTex format.

                :return rate_string0: The formatted rate constant in LaTeX form.
            """

            # The rate constant list must be a string.
            rate_constant_list0 = state.split(".")

            # Join all in uppercase letters.
            rate_string0 = "".join(rate_constant_list0)
            rate_string0 = rate_string0.upper()

            return rate_string0

        def single_term(term: tuple) -> str:
            """ Gets the string representation of a single term.

                :param term: A 2-tuple, where the first index represents the
                particle and the second index represents the site where the
                particle is sitting.

                :return: A string representation of the term.
            """

            # Format the term properly.
            term_string = f"{term[0]}{term[1]}"

            return term_string

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Auxiliary variables.
        state_string = ""

        # If it is the rate constant the one that is to be formatted.
        if rate_constant:
            return get_rate_constant()

        if order == 0 or len(state) <= order:
            # Format the string.
            state_string += get_enclosed_state(state)

        else:
            # Format the string.
            state_list = get_nth_order()
            state_string += get_nth_order_string()

        # Format the state in derivative form.
        if derivative:
            state_string = "D[" + state_string + ", t]"

        return state_string
