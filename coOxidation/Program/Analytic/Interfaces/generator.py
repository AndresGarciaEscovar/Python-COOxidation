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

        - self.order: The order to which the equations must be approximated.

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
    def order(self) -> int:
        """ Gets the order to which the equations must be approximated.

            :return: Returns the order to which the equations must be
             approximated.
        """
        return self.__order

    @order.setter
    def order(self, order: int) -> None:
        """ Sets the order to which the equations must be approximated.

            :param order: The order to which the equations must be approximated.
        """
        self.__order = min(int(abs(order)), self.sites_number)

    @order.deleter
    def order(self) -> None:
        """ Cannot delete this variable."""
        raise AttributeError("Cannot delete the order.")

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

        self.__sites_number = int(sites_number)
        self.validate_sites_number()

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
    # Methods.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Get Methods.
    # --------------------------------------------------------------------------

    @abstractmethod
    def get_associated_operations(self) -> tuple:
        """ Returns a tuple with the 3-tuples that contain the operations that
            can be applied to the states of the system, the order of the process
            and the string representation of the state.

            :return: A tuple with the 3-tuples that contain the order of the
             process, the string representation of rate associated with the
             process and the operations that can be applied to the states. Given
             in the order: (process order, process rate constant, pointer to
             function).
        """

        # Define the operations.
        process_information = (
            # (process order, process rate constant, process function)
            (None, None, None),
        )

        return process_information

    @abstractmethod
    def get_constraints(self, states: list) -> tuple:
        """ Given a set of numbered states, it gets the constraints of the
            system, i.e., the probability identities, 1 = sum(x in X) P(x),
            P(x) = sum(y in Y) P(x,y), P(y) = sum(x in X) P(x,y), etc; with
            0 <= P(x) <= 1 for x in X.

            :param states: A tuple that containst all of the "left-hand" states
             of the system.

            :return: The constrainsts of the system as equalities.
        """
        pass

    def get_contracted_state(self, states: list, index: int = -1) -> tuple:
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

            for i0, state0 in enumerate(states0):
                particles0_, indexes0_ = self.get_state_elements(state0)

                indexes0.append(indexes0_)
                particles0.append(list(particles0_))

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        # Make the index positive.
        while index < 0:
            index += len(states[0])

        # Set the variables.
        indexes, particles = [], []
        get_particles_and_indexes(states, indexes, particles)
        particles = np.array(particles)

        # Compare indexes.
        indexes = set(indexes)
        if not len(indexes) == 1:
            return tuple(), states
        indexes = list(indexes.pop())

        # Particles must match the states.
        if not set(tuple(particles[:, index])) == set(self.states):
            return tuple(), states

        # Delete the column and validate.
        particles = np.delete(particles, index, axis=1)
        particles = {tuple(particle) for particle in particles}
        if not len(particles) == 1:
            return tuple(), states
        particles = list(particles.pop())

        # Remove the given index.
        indexes.pop(index)
        if len(indexes) == 0:
            return (1,), states

        contracted = tuple(item for item in zip(particles, indexes))
        return contracted, states

    def get_decay_states(self, state: tuple, operations: dict) -> tuple:
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

        self.validate_sites_number()
        decay_states = {}
        keys = operations.keys()

        # ----------------------------------------------------------------------
        # Get the decay states.
        # ----------------------------------------------------------------------

        for key in keys:
            decay_states[key] = operations[key](state)

        return state, decay_states

    def get_is_substate(self, state1: tuple, state2: tuple) -> bool:
        """ Determines if state 1 is a substate, or proper substate, of state 2;
            they must be in the same format; order matters in this case.

            :param state1: The state that is to be found within state2.

            :param state2: The state where state1 is going to be searched.

            :return: True if state1 is a substate, or proper substate, of
             state2. False, otherwise.
        """

        if len(state1) > len(state2):
            return False

        if state1 == state2:
            return True

        substate1 = sorted(state1, key=lambda x: (x[1], x[0]))
        substate2 = list(map(tuple, itertools.combinations(state2, len(substate1))))

        for substate in substate2:
            if substate1 == sorted(substate, key=lambda x: (x[1], x[0])):
                return True

        self.validate_sites_number()

        return False

    def get_multiplicity(self, processes: dict) -> dict:
        """ Given the decay/created states dictionary, it returns a dictionary
            with the UNIQUE states for each process and their multiplicity.

            :param processes: The dictionary of processes associated with a
             state.

            :return: A dictionary with the UNIQUE states for each process and
             their multiplicity.
        """

        self.validate_sites_number()
        multiplicity = {}
        keys = processes.keys()

        for key in keys:
            unique = set(processes[key])
            multiplicity[key] = [(state, processes[key].count(state)) for state in unique]

        return multiplicity

    def get_nth_order_equations(self, display: bool = False) -> None:
        """ Gets the nth order equations for the system in the format: (state,
            decay_states, create_states); where "decay_states" and
            "create_states" are the dictionaries that contain the states to
            which the "state" decays, or the states from where the state is
            created, due to the different process in the system.

            Saves the equations to the variables for them to be processed later.

            :param display: If the table of equations must be displayed in the
             console.
        """

        # ALWAYS empty the equations list.
        self.equations = []

        # Get the states to perform the calculation.
        states_left = self.get_states_left()
        states_right = self.get_states_right()

        # Get the decay states.
        function = self.get_decay_states
        processes = self.get_process_functions()
        decay_states = list(map(lambda x: function(x, processes), states_right))

        # Get the create and decay state(s) of the left-hand states.
        for state_left in states_left:
            dictionary_create = self.get_products_create(state_left, decay_states)
            dictionary_decay = self.get_products_decay(state_left, decay_states)
            self.equations.append((state_left, dictionary_create, dictionary_decay,))

        self.get_constraints(states_left)

        if display:
            self.print_equation_states()

    @abstractmethod
    def get_numbering(self, state: tuple) -> list:
        """ Returns a tuple with the possible numbering a state of length N can
            have. The format of a SINGLE state must be given in the format:
            ( (particle_at_site1, numbering_scheme1),..., (particle_at_siteN,
            numbering_schemeN))

            :param state: The state to be numbered.

            :return: The list of possible numbered states in the given format.
        """
        pass

    def get_process_functions(self) -> dict:
        """ Returns all the pointers to the functions that operate on the
            different states to, potentially, modify them.

            :return: A dictionary, whose keys are the string represenation of
             the rates, with the functions that will potentially modify a state.
        """

        operations = self.get_associated_operations()
        process_functions = dict((operation[1], operation[2],) for operation in operations)

        return process_functions

    def get_process_orders(self) -> dict:
        """ Returns all the dictionary of integers that represent the minimum
            number of sites required for the given processes to take place.

            :return: A dictionary, whose keys are the string represenation of
             the rates, with the integers that represent the minimum number of
             sites required for the given processes to take place.
        """

        operations = self.get_associated_operations()
        rates_orders = dict((operation[1], operation[0], ) for operation in operations)

        return rates_orders

    def get_process_rates(self) -> tuple:
        """ Returns all the string representation of the rates associated with
            the class.

            :return: A tuple with the string rates associated with
             the class.
        """

        operations = self.get_associated_operations()
        rates_strings = tuple(operation[1] for operation in operations)

        return rates_strings

    def get_products_create(self, state: tuple, states_decay: list):
        """ Given a state and the list of states, whose decay states due to the
            different process must be included, it returns of the states that
            create the given state.

            :param state: The state that must be created from the other states.

            :param states_decay: A list of states and their decay products. It
             must be in the format. (state, 2-tuple with state and dictionary
             with decay states).

            :return: A 2-tuple with the state and a dictionary with the UNIQUE
             lowest order states that will create the state through a given
             process.
        """

        keys = states_decay[0][1].keys()
        create_dictionary = {key: [] for key in keys}

        for state_decay in states_decay:
            if self.get_is_substate(state, state_decay[0]):
                continue

            for key in keys:
                for state_0 in state_decay[1][key]:
                    if not self.get_is_substate(state, state_0):
                        continue

                    create_dictionary[key].append(state_decay[0])

        create_dictionary = self.reduce_to_unique_states(create_dictionary, state)
        create_dictionary = self.get_multiplicity(create_dictionary)

        return create_dictionary

    def get_products_decay(self, state: tuple, states_decay: list):
        """ Given a state and the list of states, whose decay states due to the
            different process must be included, it returns of the states that
            make the first state decay.

            :param state: The state whose decay process are to be obtained.

            :param states_decay: A list of states and their decay products. It
             must be in the format. (state, 2-tuple with state and dictionary
             with decay states).

            :return: A 2-tuple with the state and a dictionary with the UNIQUE
             lowest order states that will make it decay through a given
             process.
        """

        keys = states_decay[0][1].keys()
        decay_dictionary = {key: [] for key in keys}

        for state_decay in states_decay:
            if not self.get_is_substate(state, state_decay[0]):
                continue

            for key in keys:
                for state_0 in state_decay[1][key]:
                    if self.get_is_substate(state, state_0):
                        continue

                    decay_dictionary[key].append(state_decay[0])

        decay_dictionary = self.reduce_to_unique_states(decay_dictionary, state)
        decay_dictionary = self.get_multiplicity(decay_dictionary)

        return decay_dictionary

    def get_state_elements(self, state: tuple) -> tuple:
        """ Gets two tuples one with the elements of the given state; one with
            all the particles and the other one with the indexes.

            :param state: The representation of the state, that must be in
            the standard form.

            :return: A 2-tuple that contains the particles in the state and its
            numerical indexes in the order:

            (particles, numerical indexes)
        """

        particle_list = tuple(entry[0] for entry in state)
        index_list = tuple(entry[1] for entry in state)

        self.validate_sites_number()

        return particle_list, index_list

    def get_states(self, order: int = 1) -> list:
        """ Given the order, it returns a list of ALL the possible combinations
            of the system variables, i.e., all the possible combinations of the
            variables in N slots, where N=order; NON-NUMBERED.

            :param order: The order of the requested states.

            :return: A list of all the possible states of the given order.
        """
        return list(product(*[cp.deepcopy(self.states) for _ in range(order)]))

    def get_states_left(self) -> list:
        """ Returns a list of the numbered lowest order states for the
            equations to be written in; i.e., the states that have the
            differential operator d/dt.

            :return: A list of the non-numbered lowest order states for the
             equations to be written in.
        """

        # Setup the variables.
        states = []
        states_0 = []
        order0 = 1 if self.order == 0 else self.order

        if order0 < self.sites_number:
            for i in range(1, order0 + 1):
                states_0.extend(self.get_states(i))

        else:
            states_0.extend(self.get_states(self.sites_number))

        # For every un-numbered state.
        for state_0 in states_0:
            states.extend(self.get_numbering(state_0))

        return states

    def get_states_right(self) -> list:
        """ Returns a list of the numbered states that will potentially appear
            in the derivative term.

            :return: A list of the non-numbered lowest order states for the
             equations to be written in.
        """

        # Set the variables.
        orders = []
        states = []
        states_0 = []
        order0 = 1 if self.order == 0 else self.order

        if order0 < self.sites_number:
            processes_orders = list(set(order_0[0] for order_0 in self.get_associated_operations()))

            for order_0 in range(1, order0 + 1):
                orders.extend([order_0 + order_1 - 1 for order_1 in processes_orders])

            # Order vectors cannot be longer than the total number of sites.
            maximum_order = max(orders)
            maximum_order = min(maximum_order, self.sites_number)

            for i in range(1, maximum_order + 1):
                states_0.extend(self.get_states(i))

        else:
            states_0.extend(self.get_states(self.sites_number))

        # For every un-numbered state.
        for state0 in states_0:
            states.extend(self.get_numbering(state0))

        return states

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

            if state0 is None or state0 == tuple():
                return ""

            state0_ = state0[0]
            multiplicity0 = state0[1]

            state0_ = "".join(["<", ",".join(map(lambda x: f"{x[0]}{x[1]}", state0_)), ">"])
            state0_ = "".join([str(multiplicity0), state0_]) if multiplicity0 > 1 else state0_

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

            keys0 = sorted(column_names0.keys(), reverse=True)
            lengths0 = list(map(lambda x: len(column_names0[x]), keys0))

            for key0 in equation0[1].keys():
                states0 = ", ".join([format_state(state0) for state0 in equation0[2][key0]])
                states0_ = ", ".join([format_state(state0) for state0 in equation0[1][key0]])
                lengths0 = [max(item0) for item0 in zip(lengths0, [len(key0), len(states0), len(states0_)])]

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

            states0 = ", ".join(list(map(lambda x: format_state(x), equation0[2][key0])))
            states0_ = ", ".join(list(map(lambda x: format_state(x), equation0[1][key0])))
            strings0 = [key0, states0, states0_]

            for i0, string0 in enumerate(strings0):
                if i0 == 0:
                    strings0[i0] = f"{strings0[i0]:^{column_widths0[i0]}}"
                    continue

                strings0[i0] = f"{strings0[i0]:<{column_widths0[i0]}}"

            strings0 = "".join([get_separator(column_widths), "\n", "|", "|".join(strings0), "|", "\n"])

            return strings0

        def get_header(column_names0: dict, column_widths0: tuple) -> str:
            """ Returns the formatted header of the table.

                :param column_names0: The dictionary with the names of the
                 columns.

                :param column_widths0:The widths of each column to properly
                 print the lines and separators.

                :return: The string that represents the header of the table.
            """

            header0 = sorted(column_names0.keys(), reverse=True)
            header0 = [f"{column_names0[key0]:^{column_widths0[i0]}}" for i0, key0 in enumerate(header0)]
            header0 = "".join([get_separator(column_widths), "\n", "|", "|".join(header0), "|"])

            return header0

        def get_separator(column_widths0: tuple) -> str:
            """ Gets a string with the row separator.

                :param column_widths0: The widths of each column to properly
                 print the lines and separators.

                :return: A string with the row separator.
            """

            separator0 = ["-" * length0 for length0 in column_widths0]
            separator0 = "".join(["-", "-".join(separator0), "-"])

            return separator0

        def validate_equations() -> bool:
            """ Validates the equations list is not empty.

                :return: True, if there are equations to display. False,
                 otherwise.
            """

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

        # The dictionary with the names of the columns.
        column_names = {
            "create": "Create Processes",
            "decay": "Decay Processes",
            "key": "Rate Constant"
        }

        # Get the column contents.
        for i, equation in enumerate(self.equations):
            print("State: ", format_state((equation[0], 1)))

            column_widths = get_column_widths(column_names, equation)
            str_ = get_header(column_names, column_widths) + "\n"
            for j, key in enumerate(equation[1].keys()):
                str_ += get_equation(key, equation, column_widths)

            str_ += get_separator(column_widths) + "\n"
            print(str_)
            del str_

    # --------------------------------------------------------------------------
    # Reduce Methods.
    # --------------------------------------------------------------------------

    @abstractmethod
    def reduce_to_unique_states(self, state_list, target_state):
        """ Given a list of states, it attempts to contract

            :param state_list: The list of states to be reduced.

            :param target_state: The state that is being targeted to appear in
            the reduced list.

            :return: A list of the reduced states in the format (state,
            multiplicity).
        """
        pass

    # --------------------------------------------------------------------------
    # Save Methods.
    # --------------------------------------------------------------------------

    @abstractmethod
    def save_equations(self, file: str = None, format_type: str = None, order: int = 0, save_path: str = None) -> None:
        """ Generates the equations in the requested format, with the terms
            approximated to the given order.

            :param file: The name of the file where the equations are to be
             saved; must be extensionless. Named equations by default.

            :param format_type: A string that represents the format of the
             requested equations. NOT case sensitive, e.g., "A" = "a".

            :param order: An integer that represents the order to which the
             equations must be approximated. Zeroth order, or less, means that
             the equations will be written in an exact way.

            :param save_path: The path of the file where the equations must be
             saved.

            :param: The path where a file with the equations is to be created,
             if at all. None, by default.
        """
        pass

    # --------------------------------------------------------------------------
    # Validate Methods.
    # --------------------------------------------------------------------------

    def validate_sites_number(self) -> None:
        """ Validates that the constrainst of the system have the proper format.
        """

        if not isinstance(self.sites_number, (int,)) or self.sites_number < 1:
            raise ValueError(
                "The number of sites must an integer greater than zero."
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

        # sites_number must always go first.
        self.sites_number = sites_number
        self.states = states
        self.constraints = []
        self.equations = []
        self.order = 0
