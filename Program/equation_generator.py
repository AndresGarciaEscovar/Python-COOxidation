""" The class that serves as the base to create an equation generator.
"""

# Imports: General.
import copy as cp
import itertools
import numpy as np

from abc import ABC, abstractmethod
from collections.abc import Iterable
from itertools import product


class EquationGenerator(ABC):
    """ Generates the differential equations in different formats.

        :param self.constraint_equations: The list where the constraint
        equations will be saved.

        :param self.equations: The list where the equations will be saved.

        :param self.sites: The maximum number of sites.

        :param  self.states: The UNIQUE states in which each site can be in.
    """

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Getters, Setters and Deleters.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    @property
    def constraint_equations(self):
        """ Gets the equations of the system.
        """
        return self.__constraint_equations

    @constraint_equations.setter
    def constraint_equations(self, constraint_equations):
        """ Sets the number of sites in the lattice.
        """

        # Verify that the parameter is a list.
        if not isinstance(constraint_equations, (list,)):
            raise ValueError("The constraint equations variable must be a list.")

        self.__constraint_equations = constraint_equations

    @constraint_equations.deleter
    def constraint_equations(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def equations(self):
        """ Gets the equations of the system.
        """
        return self.__equations

    @equations.setter
    def equations(self, equations):
        """ Sets the number of sites in the lattice.
        """

        # Verify that the parameter is a list.
        if not isinstance(equations, (list,)):
            raise ValueError("The equations variable must be a list.")

        self.__equations = equations

    @equations.deleter
    def equations(self):
        pass

    # --------------------------------------------------------------------------

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
    # Get Methods.
    # --------------------------------------------------------------------------

    @abstractmethod
    def save_equations(self, file_name="equations", format_type="latex", order=0, save_path=None):
        """ Generates the equations in the requested format, with the terms
            approximated to the given order.

            :param file_name: The name of the file where the equations are to be
            saved; must be extensionless. Named equations by default.

            :param format_type: A string that represents the format of the requested
            equations. NOT case sensitive, e.g., "A" = "a".

            :param order: The order to which the equations must be approximated.
            Zeroth order, or less, means that the equations will be written in
            an exact way.

            :param save_path: The path where a file with the equations is to be
            created, if at all. None, by default.

            :param together: If the equations must be printed in a stand-alone
            format of gathered format.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------
        pass

    def get_nth_order_equations(self, order=0, print_equations=False):
        """ Gets the nth order equations for the system in the format
                (state, decay_states, create_states),
            where "decay_states" and "create_states" are the dictionaries that
            contain the states to which the "state" decays, or the states from
            where the state is created, due to the different process in the
            system.

            :param order: The lowest order to which the equations must be
            given.

            :param print_equations: If the table of equations must be printed.
        """

        # ----------------------------------------------------------------------
        # ALWAYS
        # ----------------------------------------------------------------------

        # Empty the equations list.
        self.equations = []

        # ----------------------------------------------------------------------
        # Get the states to perform the calculation.
        # ----------------------------------------------------------------------

        # Get the lowest order states.
        states_left_hand = self._get_states_left(order)

        # Get the other involved states.
        states_right_hand = self._get_states_right(order)

        # ----------------------------------------------------------------------
        # Get the decay states.
        # ----------------------------------------------------------------------

        # An alias to the function to make it shorter.
        decay_func = self._get_decay_states

        # Dictionary of processes.
        process_dict = self._get_process_functions()

        # Get the decay states for each of the right-hand states.
        decay_states = list(map(lambda x: decay_func(x, process_dict), states_right_hand))

        # ----------------------------------------------------------------------
        # Get the decay state(s) of the left-hand states.
        # ----------------------------------------------------------------------

        # Get both the decay and creation dictionaries for each state.
        for state_left_hand in states_left_hand:
            # Get the list of decay states.
            decay_dictionary = self._get_products_decay(state_left_hand, decay_states)

            # Get the list of create states.
            create_dictionary = self._get_products_create(state_left_hand, decay_states)

            # Append it to the equations.
            self.equations.append((state_left_hand, create_dictionary, decay_dictionary,))

        # Get the constraints.
        self._get_constraint_equations(states_left_hand)

        # Print the equations if needed.
        _ = self.print_equation_states() if print_equations else None

    # --------------------------------------------------------------------------
    # Other Methods.
    # --------------------------------------------------------------------------

    def print_equation_states(self):
        """ Prints the states of the equations to the screen.
        """

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Check that the equations list is not empty.
        if len(self.equations) == 0:
            print("There are no equations to print.")

        # Get the keys.
        keys = [str(key) for key in self._get_process_functions().keys()]

        # Auxiliary variables.
        string_create = "Create Processes"
        string_decay = "Decay Processes"
        string_key = "Key"

        # For every equation in the equations list.
        for equation in self.equations:

            # Print the state.
            print("State: ", equation[0])

            # Get the width of the key column.
            cw1 = max(len(string_key), max(map(len, keys)))

            # Get the width of the create process column.
            tmp_create = {key: ", ".join([str(state) for state in equation[1][key]]) for key in keys}
            cw2 = max(len(string_create), max(map(len, [tmp_create[key] for key in keys])))

            # Get the width of the decay process column.
            tmp_decay = {key: ", ".join([str(state) for state in equation[2][key]]) for key in keys}
            cw3 = max(len(string_decay), max(map(len, [tmp_decay[key] for key in keys])))

            # Format the separator string.
            decay_string = ("|", ("-" * cw1), "|", ("-" * cw2), "|", ("-" * cw3), "|",)

            # Print the header.
            print(*decay_string)
            print("|", f"{string_key:^{cw1}}", "|", f"{string_create:^{cw2}}", "|", f"{string_decay:^{cw3}}", "|", )
            print(*decay_string)

            # For every key.
            for key in keys:
                # Print the table entry.
                print("|", f"{key:^{cw1}}", "|", f"{tmp_create[key]:<{cw2}}", "|", f"{tmp_decay[key]:<{cw3}}", "|", )

                # Remember to print the separator.
                print(*decay_string)

            print("")

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Private Interface.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Generate Methods.
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # Get Methods.
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

    @abstractmethod
    def _get_constraint_equations(self, states):
        """ Given a set of numbered states, it gets the constraints of the
            system, i.e., the probability identities, 1 = sum(x in X) P(x),
            P(x) = sum(y in Y) P(x,y), P(y) = sum(x in X) P(x,y), etc; with
            0 <= P(x) <= 1 for x in X.

            :param states: ALL of the "left-hand" states of the system.

            :return:  The constrainst of the system as equalities.
        """
        pass

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

        # Get the states up to the given order.
        if order0 < self.number_of_sites:
            # Up until the requested order.
            for i in range(1, order0 + 1):
                # Get the un-numbered states.
                states_0.extend(self._get_states(i))
        else:
            # Exact Equations.
            states_0.extend(self._get_states(self.number_of_sites))

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
        states_0 = []

        # Fix the order, if needed.
        order0 = 1 if order == 0 else order
        order0 = order0 if order < self.number_of_sites else self.number_of_sites

        # Get the states up to the given order.
        if order0 < self.number_of_sites:
            # Get the orders of the process.
            processes_orders = list(set(order_0[0] for order_0 in self._get_associated_operations()))

            # For each length of states.
            for order_0 in range(1, order0 + 1):
                # Get the maximum order.
                orders.extend([order_0 + order_1 - 1 for order_1 in processes_orders])

            # Order vectors cannot be longer than the total number of sites.
            maximum_order = max(orders)
            maximum_order = min(maximum_order, self.number_of_sites)

            # For all the orders, up to the maximum order.
            for i in range(1, maximum_order + 1):
                # Get the states.
                states_0.extend(self._get_states(i))

        else:
            # Exact Equations.
            states_0.extend(self._get_states(self.number_of_sites))

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

        # List where the constraint equations are saved.
        self.constraint_equations = []

        # List where the equations are saved.
        self.equations = []

    # --------------------------------------------------------------------------
    # Dunder Methods.
    # --------------------------------------------------------------------------
