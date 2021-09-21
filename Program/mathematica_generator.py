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

    def generate_equations(self, gather_by_state=False, format_string="latex", order=0, save_file_name="equations"):
        """ Generates the equations in LaTeX format.

            :param gather_by_state: True, if the terms must be factorized by
            state. False, if the terms must be factorized by reaction constant.

            :param format_string: The type of formatting that must be applied to
            the strings. LaTeX is the default setting.

            :param order: The order to which approximate the terms. Must
            be an integer greater than or equal to zero.

            :param save_file_name: The name of the file where the strings will
            be saved.

            :return: A string that represents the equations in LaTeX format.
        """

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # List where the output strings are to be saved.
        output_strings = []

        # ----------------------------------------------------------------------
        # Get the equations.
        # ----------------------------------------------------------------------

        # Get all the strings.
        for equation in self.equations:
            # Save the output string to the list to be saved, if needed.
            output_strings.append(self._generate_equation(equation, gather_by_state, format_string, order))

        # ----------------------------------------------------------------------
        # Save the strings to a file.
        # ----------------------------------------------------------------------

        # Join the equations to save.
        to_save = "\n\n".join(output_strings)

        # Save it to the requested file.
        with open(save_file_name, "w") as fl:
            fl.write(to_save)

    # --------------------------------------------------------------------------
    # Get methods.
    # --------------------------------------------------------------------------

    def get_0th_order_equations(self, print_equations=False):
        """ Gets the lowest order exact equations.

            :param print_equations: True, if the equations are to be printed
            once the process finishes. False, otherwise.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def empty_dictionaries():
            """ Resets the dictionaries to empty dictionaries.
            """

            # Set all the entries in the dictionary to empty.
            for key0 in keys:
                decay_dictionary[key0] = []
                create_dictionary[key0] = []

        def get_creation_states(state0, state1):
            """ Appends the states to which the state0 decays for each process.

                :param state0: The list that contains the information of the
                state0, that is the state which is tested for decay.

                [original_state, decay_dictionary, create_dictionary]

                :param state1: The list that contains the information of the
                state which must contain state0.

                [original_state, decay_dictionary]
            """

            # Check all the states belonging to the specific key.
            for key0 in keys:
                # Check all the states generated by the process.
                for statei in state1[1][key0]:
                    # Check that the state is a substate of the state, in the lowest order possible.
                    if len(statei) == process_orders[key0] and self._get_state1_in_state2(state0[0], statei):
                        # Add the state to the list of generated states.
                        state0[2][key0].append(state1[0])

        def get_decay_states(state0, state1):
            """ Appends the states to which the state0 decays for each process.

                :param state0: The list that contains the information of the
                state0, that is the state which is tested for decay.

                [original_state, decay_dictionary, create_dictionary]

                :param state1: The list that contains the information of the
                state which must contain state0.

                [original_state, decay_dictionary]
            """

            # Only proceed if the state is in the process.
            if not self._get_state1_in_state2(state0[0], state1[0]):
                return

            # Check all the keys.
            for key1 in keys:
                # Append the states to the dictionary if needed.
                if len(state1[0]) == process_orders[key1] and len(state1[1][key1]) > 0:
                    state0[1][key1].extend([cp.deepcopy(state1[0]) for _ in state1[1][key1]])

        def get_involved_states():
            """ Gets the states that are involved in the calculation.

                :return involved_states0: A list with the states that are
                involved in the calculation.
            """

            # Auxiliary variables.
            involved_states0 = []

            # The states of the orders that show in the lowest order equations.
            for order0 in self._get_orders():
                involved_states0.extend(self._get_states(order0))

            return involved_states0

        def get_resulting_states():
            """ Gets the resultant states after ALL processes are applied to ALL
                the states.
            """

            # Auxiliary variables.
            resultant_states0 = []

            for state0 in involved_states:
                # Empty the dictionary first.
                decay_dictionary_state0 = {}

                # Get the states to which each state decays due to a process.
                for i, process0 in enumerate(process_functions):
                    decay_dictionary_state0[keys[i]] = process0[2](state0)

                # Append the results.
                resultant_states0.append([state0, decay_dictionary_state0])

            return resultant_states0

        def print_states():
            """ Prints a table of the decay states and the creation states of
                all the states in the system.
            """

            # The titles of the columns.
            keyis = "Process"
            decay = "Decay States"
            create = "Create States"

            # Get the information on the strings.
            tmp_key = [len(key) for key in keys]
            tmp_state_create = [len(create)]
            tmp_state_decay = [len(decay)]
            for state0 in self.equations:
                for key in keys:
                    tmp_state_decay.append(len(str(state0[1][key])))
                    tmp_state_create.append(len(str(state0[2][key])))
            tmp_key.append(len(keyis))

            # Get the maximum lengths.
            max_create = max(tmp_state_create)
            max_decay = max(tmp_state_decay)
            max_key = max(tmp_key)

            # Format strings.
            str1 = "".join(["-" for _ in range(max_key)])
            str2 = "".join(["-" for _ in range(max_decay)])
            str3 = "".join(["-" for _ in range(max_create)])

            # Print the formatted states.
            for state0 in self.equations:
                print(f"Initial state is {state0[0]}:")
                print(f"\t", "|", str1, "|", str2, "|", str3, "|")
                print(f"\t", "|", f"{keyis:^{max_key}}", "|", f"{decay:^{max_decay}}", "|", f"{create:^{max_create}}", "|")
                print(f"\t", "|", str1, "|", str2, "|", str3, "|")
                for key0 in keys:
                    print(f"\t", "|", f"{key0:<{max_key}}", "|", f"{str(state0[1][key0]):<{max_decay}}", "|", f"{str(state0[2][key0]):<{max_create}}", "|")
                    print(f"\t", "|", str1, "|", str2, "|", str3, "|")
                print("")

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # First, empty the equations list.
        self.equations = []

        # Get the information of the process functions and process orders.
        process_functions, process_orders = self._get_process_functions()

        # Get the lowest state variables, with the numbering.
        lowest_states = self._get_states(1)

        # Get ALL the states that are potentially involved in the calculation.
        involved_states = get_involved_states()

        # ----------------------------------------------------------------------
        # Get the generated states by each process.
        # ----------------------------------------------------------------------

        # Auxiliary variables.
        decay_dictionary = {}
        create_dictionary = {}

        # Get the keys.
        keys = self._get_keys()

        # Get the resulting states from making the operations.
        resultant_states = get_resulting_states()

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Empty the decay and creation dictionaries.
        empty_dictionaries()

        # Check how each state can be generated.
        for low_state in lowest_states:
            # Set the format of the state to the final state.
            low_state = [low_state, cp.deepcopy(decay_dictionary), cp.deepcopy(create_dictionary)]

            # Check every state.
            for resultant_state in resultant_states:
                # Get the decay states.
                get_decay_states(low_state, resultant_state)

                # If the states are different we testing creation.
                get_creation_states(low_state, resultant_state)

            # Add the equation to the list.
            self.equations.append(cp.deepcopy(low_state))

        # Print the table with all the states.
        if print_equations:
            print_states()

    def get_nth_order_equations(self, order=0, print_equations=False):
        """ Gets the equations up to the nth order term, i.e., the equations
            of states from 0 to n-1 order. If the order is the same as the
            number of sites, it gives the EXACT equations for the system.

            :param order: The order to which the equations must be obtained. It
            must be an integer in the range [0, self.number_of_sites].

            :param print_equations: True, if the equations are to be printed
            once the process finishes. False, otherwise.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def empty_dictionaries():
            """ Resets the dictionaries to empty dictionaries.
            """

            # Set all the entries in the dictionary to empty.
            for key0 in keys:
                decay_dictionary[key0] = []
                create_dictionary[key0] = []

        def get_creation_states(state0, state1):
            """ Appends the states that create state0 for each process.

                :param state0: The list that contains the information of the
                state0, that is the state which is tested for creation.

                :param state1: The list that contains the information of the
                state which must contain state0.

                [original_state, decay_dictionary]
            """

            # A state cannot be created by itself.
            if self._get_state1_in_state2(state0, state1[0]):
                print("Rejected:")
                print(state0)
                print(state1[0])
                print("")
                return

            print("Accepted:")
            print(state0)
            print(state1[0])
            print("")

            # # Go through all the process the keys.
            # for key0 in keys:
            #     # Go through all the creation states due to a process.
            #     for state0_0 in state1[1][key0]:
            #         # Only append if the original state has decayed.
            #         if self._get_state1_in_state2(state0, state0_0):
            #             create_dictionary[key0].append(cp.deepcopy(state1[0]))
            #             break

        def get_decay_states(state0, state1):
            """ Appends the states to which the state0 decays for each process.

                :param state0: The list that contains the information of the
                state0, that is the state which is tested for decay.

                :param state1: The list that contains the information of the
                state which must contain state0.

                [original_state, decay_dictionary]
            """

            # Only proceed if the state is in the process.
            if not self._get_state1_in_state2(state0, state1[0]):
                return

            # Go through all the process the keys.
            for key0 in keys:
                # Go through all the decay states due to a process.
                for state0_0 in state1[1][key0]:
                    # Only append if the original state has decayed.
                    if not self._get_state1_in_state2(state0, state0_0):
                        decay_dictionary[key0].append(cp.deepcopy(state1[0]))

        def get_filtered_contracted_states(state0, states0):
            """ Gets a list of all the UNIQUE states, with the number of times
                it appears in the list.

                :param state0: The target state to be obtained by contraction.

                :param states0: The list of states to be filtered.

                :return filtered_states0: The list of filtered states.
            """

            print(state0)

            # Proceed only if the list length is greater than zero.
            if len(states0) == 0:
                return states0

            # Get a list of all the UNIQUE states and get the count of them.
            filtered_states0 = cp.deepcopy(states0)
            filtered_states0 = sorted(list(set(filtered_states0)), key=lambda x: (len(x), x[0][0], x[0][1]))

            # Number of unique states.
            unique_state_length0 = -1

            # Contract the states if needed.
            while not unique_state_length0 == len(filtered_states0):
                # Get the list of all the combinations needed to contract a state.
                nplets_combinations0 = tuple(itertools.combinations(filtered_states0, len(self.states)))

                # If there are no triplet combinations.
                if len(nplets_combinations0) == 0:
                    break

                # Try contract the nplets at both sides.
                for nplets_combination0 in nplets_combinations0:
                    contracted0_0 = self._get_contracted_state(nplets_combination0, -1)
                    contracted0_1 = self._get_contracted_state(nplets_combination0, 0)

                    # If there are contracted indexes towards the right.
                    if len(contracted0_0[0]) > 1 and self._get_state1_in_state2(state0, contracted0_0[0]):
                        filtered_states0 = [state0_2 for state0_2 in filtered_states0 if state0_2 not in contracted0_0[1]]

                        # Add the state as needed.
                        if contracted0_0[0] not in filtered_states0:
                            filtered_states0.append(contracted0_0[0])

                        break

                    # If there are contracted indexes towards the left.
                    if len(contracted0_1[0]) > 1 and self._get_state1_in_state2(state0, contracted0_1[0]):
                        filtered_states0 = [state0_2 for state0_2 in filtered_states0 if state0_2 not in contracted0_1[1]]

                        # Add the state as needed.
                        if contracted0_1[0] not in filtered_states0:
                            filtered_states0.append(contracted0_1[0])

                        break

            # ------------------------------------------------------------------
            # Get the multiplicity of the terms.
            # ------------------------------------------------------------------

            # Auxiliary variables.
            tmp_filtered = cp.deepcopy(filtered_states0)

            # Empty the filtered states.
            filtered_states0 = []

            # Get the new filtered states.
            for state0 in tmp_filtered:
                # Get each state's multiplicity.
                num = 1 if states0.count(state0) == 0 else states0.count(state0)

                # Add to the list.
                filtered_states0.append(tuple([state0, num]))

            return filtered_states0

        def get_involved_states():
            """ Gets the states that are involved in the calculation.

                :return involved_states0: A list with the states that are
                involved in the calculation.
            """

            # Auxiliary variables.
            involved_orders0 = [order]
            involved_states0 = []
            orders0 = []

            # ------------------------------------------------------------------
            # Get the orders involved.
            # ------------------------------------------------------------------

            # Get the orders of the processes.
            orders0.extend(self._get_orders())

            # Get the order of the states for which the gains/decays will happen.
            length_states0 = set(tuple(map(len, lowest_states)))

            # Get the orders of the states involved in the calculations.
            for order0 in orders0:
                # State of the current order will always be involved.
                involved_orders0.append(order0)

                # Get the specific order of each state.
                for length_state0 in length_states0:
                    # Use the criteria to add the states.
                    involved_orders0.extend([k for k in range(order0, order0 + length_state0)])

            # Filter the involved orders.
            involved_orders0 = set(involved_orders0)
            involved_orders0 = [k for k in involved_orders0 if order <= k <= self.number_of_sites]

            # Get ALL the involved states.
            for order0 in involved_orders0:
                involved_states0.extend(self._get_states(order0))

            return involved_states0

        def get_lowest_states():
            """ Gets the states that are involved in the calculation.

                :return involved_states0: A list with the states that are involved in the
                calculation.
            """

            # Auxiliary variables.
            involved_states0 = []

            # ----------------------------------------------------------------------------------------------------------
            # Get the states that are involved in the calculation.
            # ----------------------------------------------------------------------------------------------------------

            # If the order is greater than the number of sites.
            if order >= self.number_of_sites:
                # Only get the states at the highest order possible.
                for involved_state0 in self._get_states(self.number_of_sites):
                    involved_states0.append(involved_state0)

                return involved_states0

            # Get ALL the states up to the lowest order possible.
            for k in range(2, order + 1):
                for involved_state0 in self._get_states(k):
                    involved_states0.append(involved_state0)

            return involved_states0

        def get_resulting_states():
            """ Gets the resultant states after ALL processes are applied to ALL
                the states.
            """

            # Auxiliary variables.
            resultant_states0 = []

            for state0 in involved_states:
                # Empty the dictionary first.
                decay_dictionary_state0 = {}

                # Get the states to which each state decays due to a process.
                for k, process0 in enumerate(process_functions):
                    decay_dictionary_state0[keys[k]] = process0[2](state0)

                # Append the results.
                resultant_states0.append([state0, decay_dictionary_state0])

            return resultant_states0

        def validate_order():
            """ Makes sure that the chosen order is valid, i.e., a positive
                integer greater than zero.
            """

            # The order must be a positive integer greater than zero.
            if not isinstance(order, (int,)) or order <= 0:
                raise ValueError(f"The order parameter must be an integer number greater than zero. Current value: {order}.")

        # # ----------------------------------------------------------------------
        # # Implementation.
        # # ----------------------------------------------------------------------
        #
        # # Get the zeroth order equations, these will serve as the basis.
        # if order == 0 or order < self.number_of_sites:
        #     self.get_0th_order_equations(print_equations)
        #
        #     # Get the equations to the lowest order.
        #     if order == 0 or order == 1:
        #         return

        # ----------------------------------------------------------------------
        # Get the equations for the nth order term; do not empty the
        # exact equations, these will serve as the basis.
        # ----------------------------------------------------------------------

        # Validate that the order is correct.
        validate_order()

        # Get the information of the process functions and process orders.
        process_functions, process_orders = self._get_process_functions()

        # Get states for which the equations will be obtained.
        lowest_states = get_lowest_states()

        # Get ALL the states that are potentially involved in the calculation.
        involved_states = get_involved_states()
        #
        # # ----------------------------------------------------------------------
        # # Get the generated states by each process.
        # # ----------------------------------------------------------------------

        # Auxiliary variables.
        decay_dictionary = {}
        create_dictionary = {}

        # Get the keys.
        keys = self._get_keys()

        # Get the resulting states from making the operations.
        resulting_states = get_resulting_states()

        # Go through each lowest order state.
        for i, state_0 in enumerate(lowest_states):

            # Empty the dictionaries.
            empty_dictionaries()

            # Get the decay states.
            # for j, state_1 in enumerate(resulting_states):
                # Get the decay states for the lowest order state.
                # get_decay_states(state_0, state_1)

                # Get the create states for the lowest order state.
                # get_creation_states(state_0, state_1)
                # pass

            for key in keys:
                # decay_dictionary[key] = get_filtered_contracted_states(state_0, decay_dictionary[key])
                # create_dictionary[key] = get_filtered_contracted_states(state_0, cp.deepcopy(create_dictionary[key]))
                pass

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Private Interface.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Generate methods.
    # --------------------------------------------------------------------------

    def _generate_equation(self, equation, gather_by_state, format_string="latex", order=0):
        """ Generates a single equation in the required format; ideal to
            generate LaTeX style equations for a Jupyter notebook.

            :param equation: A list that represents an equation. Must be in the
            format

            [state, decay_dictionary, creation_dictionary]

            :param gather_by_state: True, if the terms must be factorized by
            state. False, if the terms must be factorized by reaction constant.

            :param format_string: The type of formatting that must be applied to
            the strings.

            :param order: The order to which approximate the terms. Must
            be an integer greater than or equal to zero.

            :return: A string that represents the equations in LaTeX format.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def get_gathered_by_constant():
            """ Returns the LaTeX equation string gathered by rate constant.

                :return: The LaTeX equation string gathered by rate constant.
            """

            # ------------------------------------------------------------------
            # Implementation.
            # ------------------------------------------------------------------

            # Auxiliary variables.
            equation_string0 = ""

            # Go through each rate constant.
            for key0 in keys:
                # Get the representation of the key.
                key_string0 = get_string(key0, rate_constant=True)

                # --------------------------------------------------------------
                # Decay terms.
                # --------------------------------------------------------------

                # Get all the decay terms.
                decay_terms0 = []
                if len(equation[1][key0]) > 0:
                    decay_terms0.extend(list(map(get_string, equation[1][key0])))

                # --------------------------------------------------------------
                # Creation terms.
                # --------------------------------------------------------------

                # Get all the creation terms.
                creation_terms0 = []
                if len(equation[2][key0]) > 0:
                    creation_terms0.extend(list(map(get_string, equation[2][key0])))

                # --------------------------------------------------------------
                # Join the two strings; three cases to contemplate.
                # --------------------------------------------------------------

                # There are decay terms, but no creation terms.
                if len(decay_terms0) > 0 and len(creation_terms0) == 0:
                    # Join all the terms.
                    decay_string0 = f"-{key_string0} (" + "+".join(decay_terms0) + ")"

                    # Remove parentheses if needed.
                    if len(decay_terms0) == 1:
                        decay_string0 = "".join(decay_string0.split("("))
                        decay_string0 = "".join(decay_string0.split(")"))

                    # Be sure to delete leading and trailing spaces.
                    decay_string0 = decay_string0.strip()

                    # Join the terms to the string.
                    equation_string0 += decay_string0

                # There are no decay terms, but there are creation terms.
                elif len(decay_terms0) == 0 and len(creation_terms0) > 0:
                    # Join all the terms.
                    create_string0 = f"{key_string0} (" + "+".join(creation_terms0) + ")"

                    # Remove parentheses if needed.
                    if len(creation_terms0) == 1:
                        create_string0 = "".join(create_string0.split("(")).strip()
                        create_string0 = "".join(create_string0.split(")")).strip()

                    # Be sure to delete leading and trailing spaces.
                    create_string0 = create_string0.strip()

                    # Join the terms to the string.
                    equation_string0 += "+" + create_string0

                # There are both decay terms and creation terms.
                elif len(decay_terms0) > 0 and len(creation_terms0) > 0:
                    # Set the digit string to empty.
                    decay_string0 = "-".join(decay_terms0)
                    create_string0 = "+".join(creation_terms0)

                    # Join the terms to the string.
                    equation_string0 += f"+{key_string0} (" + create_string0 + "-" + decay_string0 + ")"

            # ------------------------------------------------------------------
            # Finish formatting the equation.
            # ------------------------------------------------------------------

            # Fix the equation.
            equation_string0 = equation_string0[1:] if equation_string0[0] == "+" else equation_string0
            equation_string0 = " + ".join(equation_string0.split("+")).strip()
            equation_string0 = " - ".join(equation_string0.split("-")).strip()
            equation_string0 = "-" + equation_string0[2:] if equation_string0[0] == "-" else equation_string0

            return equation_string0

        def get_gathered_by_state():
            """ Returns the LaTeX equation string gathered by state.

                :return: The LaTeX equation string gathered by state.
            """

            # ------------------------------------------------------------------
            # Auxiliary functions.
            # ------------------------------------------------------------------

            def get_all_states():
                """ Returns a list of all the unique states that appear in the
                    equation.

                    :return states1: A list of all the unique states that appear
                    in the equation.
                """

                # Auxiliary variables.
                states1 = []

                # Get all the states for all the processes.
                for key1 in keys:
                    # Get the states that make the current state decay.
                    for state1 in equation[1][key1]:
                        states1.append(state1)

                    # Get the states create the current state.
                    for state1 in equation[2][key1]:
                        states1.append(state1)

                # Get a list of unique states.
                states1 = list(set(tuple(states1)))

                return states1

            def get_key_count(state1, key1):
                """ Returns the number of overall number of times the state
                    appears for a certain rate constant.

                    :param state1: The state that is being looked for.

                    :param key1: The key of the entry in the dictionary.

                    :return: The number of overall number of times the state
                    appears for a certain rate constant.
                """

                # Count the number of times the constant appears in the
                # each of the keys.
                key_count1 = 0

                # Subtract one if the state decays.
                if state1 in equation[1][key1]:
                    key_count1 -= 1

                # Subtract one if the state is created.
                if state1 in equation[2][key1]:
                    key_count1 += 1

                return key_count1

            # ------------------------------------------------------------------
            # Implementations.
            # ------------------------------------------------------------------

            # Auxiliary variables.
            equation_string0 = ""

            # Get all the unique states..
            states0 = get_all_states()

            # Look for the state in the original dictionaries.
            for state0 in states0:
                # Save the keys that are present.
                keys_present0 = []

                # Look for the state in all the keys.
                for key0 in keys:
                    # Get the formatted key.
                    key_string0 = get_string(key0, rate_constant=True)

                    # Get the number of times the state appears in the terms.
                    key_count0 = get_key_count(state0, key0)

                    # If the state does not exists in either of the given key,
                    # go the next one.
                    if key_count0 == 0:
                        continue

                    # Format the key properly.
                    key_string0 = str(abs(key_count0)) + key_string0 if abs(key_count0) > 1 else key_string0
                    key_string0 = "-" + key_string0 if key_count0 < 0 else key_string0

                    # Add the key to the list of keys
                    keys_present0.append(key_string0)

                # --------------------------------------------------------------
                # Join to the equation string.
                # --------------------------------------------------------------

                if len(keys_present0) == 1:
                    # Format the string properly.
                    key_string0_1 = "".join(keys_present0) + " " + get_string(state0)

                    # Add to the equation string.
                    equation_string0 += key_string0_1 if key_string0_1[0] == "-" else "+" + key_string0_1
                else:
                    # Format the string properly.
                    key_string0_1 = "".join([string0 if string0[0] == "-" else "+" + string0 for string0 in keys_present0]).strip()
                    key_string0_1 = key_string0_1[1:] if key_string0_1[0] == "+" else key_string0_1
                    key_string0_1 = "(" + key_string0_1 + ") " + get_string(state0)

                    # Add to the equation string.
                    equation_string0 += "+" + key_string0_1

            # Format the equation further.
            equation_string0 = equation_string0.strip()
            equation_string0 = equation_string0[1:] if equation_string0[0] == "+" else equation_string0
            equation_string0 = " + ".join(equation_string0.split("+")).strip()
            equation_string0 = " - ".join(equation_string0.split("-")).strip()
            equation_string0 = "-" + equation_string0[2:] if equation_string0[0] == "-" else equation_string0

            return equation_string0

        def get_string(state0, rate_constant=False):
            """ Gets the string representation of the given state, or rate
                constant.

                :param state0: The state for which the string representation is
                 required.

                :param rate_constant: If it is a rate constant the one that is
                required.

                :return: The string representation of the given state in the
                requested format.
            """

            return formatter(state0, order=order, rate_constant=rate_constant)

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Define the format in which the equation must be obtained.
        formatter = EquationFormatter.get_formatter(format_string)

        # Get the keys of the processes.
        keys = self._get_keys()

        # Get the differential of the state.
        equation_string = formatter(equation[0], derivative=True, order=order) + " ="
        equation_string += "= " if format_string == "mathematica" else " "

        # Write the results gathered in the requested order.
        equation_string += get_gathered_by_state().strip() if gather_by_state else get_gathered_by_constant().strip()

        return equation_string

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

            :return new_dictionary: A dictionary with the UNIQUE states for each
            process and their multiplicity.
        """

        # Auxiliary variables.
        keys = state_dictionary.keys()

        # Define the dictionary that contains the UNIQUE processes and their
        # multiplicity.
        new_dictionary = {}

        print(state_dictionary)

        # For each process.
        for key in keys:
            # Get the unique states.
            unique_states = set(state_dictionary[key])


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

            :param sites: The number of sites that the system has.

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
