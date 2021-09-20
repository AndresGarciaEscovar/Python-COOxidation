""" Writes the equations for carbon monoxide - oxygen associative reaction on
    ruthenium (111); J. Chem. Phys. 143, 204702 (2015).
    https://doi.org/10.1063/1.4936354
"""

# Imports.
import copy as cp
import itertools

import numpy as np

from collections.abc import Iterable
from itertools import product


class EquationGenerator:
    """ Generates the differential equations in different formats. Currently
        only Mathematica and LaTeX are supported.

        :param self.k_coo_er: String that represents the rate constants of
        carbon monoxide - oxygen gaseous associative desorption of surface
        oxygen; i.e., Elay-Rideal reaction

        :param self.k_coo_lh: String that represents the rate constants of
        carbon monoxide - oxygen neighboring pair associative desorption;
        i.e., Langmuir-Hinshelwoodd reaction.

        :param self.k_co_ads: String that represents the rate constants of
        adsorption of carbon monoxide.

        :param self.k_co_des: String that represents the rate constants of
        desorption of carbon monoxide.

        :param self.k_co_dif: String that represents the rate constants of
        diffusion of carbon monoxide.

        :param self.k_o_ads: String that represents the rate constants of
        adsorption of oxygen.

        :param self.k_o_des: String that represents the rate constants of
        desorption of oxygen.

        :param self.k_o_dif: String that represents the rate constants of
        diffusion of oxygen.

        :param self.o_coo_er: An integer that represents the minimum length of
        the state for carbon monoxide - oxygen gaseous associative desorption of
        surface oxygen; i.e., Elay-Rideal reaction.

        :param self.o_coo_lh: An integer that represents the minimum length of
        the state for carbon monoxide - oxygen neighboring pair associative
        desorption; i.e., Langmuir-Hinshelwoodd reaction.

        :param self.o_co_ads: An integer that represents the minimum length of
        the state for carbon monoxide adsorption.

        :param self.o_co_des: An integer that represents the minimum length of
        the state for carbon monoxide desorption.

        :param self.o_co_dif: An integer that epresents the minimum length of the state for
        carbon monoxide diffusion.

        :param self.o_o_ads: An integer that represents the minimum length of
        the state for oxygen adsoprtion.

        :param self.o_o_des: An integer tnat represents the minimum length of
        the state for oxygen desoprtion.

        :param self.o_o_dif: An integer that represents the minimum length of
        the state for oxygen diffusion.

        :param self.sites: The maximum number of sites

        :param  self.states:
    """

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Getters, Setters and Deleters.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    @property
    def k_o_ads(self):
        """ Gets the oxygen adsorption order.
        """
        return cp.deepcopy(self.__k_o_ads)

    @k_o_ads.setter
    def k_o_ads(self, k_o_ads):
        """ Sets the oxygen adsorption parameter.
        """
        self.__k_o_ads = "k.O.ads"

    @k_o_ads.deleter
    def k_o_ads(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def k_o_des(self):
        """ Gets the oxygen desorption order.
        """
        return cp.deepcopy(self.__k_o_des)

    @k_o_des.setter
    def k_o_des(self, _):
        """ Sets the oxygen desorption parameter.
        """
        self.__k_o_des = "k.O.des"

    @k_o_des.deleter
    def k_o_des(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def k_o_dif(self):
        """ Gets the oxygen diffusion order.
        """
        return cp.deepcopy(self.__k_o_dif)

    @k_o_dif.setter
    def k_o_dif(self, _):
        """ Sets the oxygen diffusion parameter.
        """
        self.__k_o_dif = "k.O.dif"

    @k_o_dif.deleter
    def k_o_dif(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def k_co_ads(self):
        """ Gets the carbon monoxide adsorption order.
        """
        return cp.deepcopy(self.__k_co_ads)

    @k_co_ads.setter
    def k_co_ads(self, _):
        """ Sets the carbon monoxide adsorption parameter.
        """
        self.__k_co_ads = "k.CO.ads"

    @k_co_ads.deleter
    def k_co_ads(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def k_co_des(self):
        """ Gets the carbon monoxide desorption order.
        """
        return cp.deepcopy(self.__k_co_des)

    @k_co_des.setter
    def k_co_des(self, _):
        """ Sets the carbon monoxide desorption parameter.
        """
        self.__k_co_des = "k.CO.des"

    @k_co_des.deleter
    def k_co_des(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def k_co_dif(self):
        """ Gets the carbon monoxide diffusion order.
        """
        return cp.deepcopy(self.__k_co_dif)

    @k_co_dif.setter
    def k_co_dif(self, _):
        """ Sets the carbon monoxide diffusion parameter.
        """
        self.__k_co_dif = "k.CO.dif"

    @k_co_dif.deleter
    def k_co_dif(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def k_coo_lh(self):
        """ Gets the carbon monoxide - oxygen reaction order.
        """
        return cp.deepcopy(self.__k_coo_lh)

    @k_coo_lh.setter
    def k_coo_lh(self, _):
        """ Sets the carbon monoxide - oxygen reaction parameter.
        """
        self.__k_coo_lh = "k.COO.lh"

    @k_coo_lh.deleter
    def k_coo_lh(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def k_coo_el(self):
        """ Gets the carbon monoxide - oxygen  gas-phase reaction order.
        """
        return cp.deepcopy(self.__k_coo_el)

    @k_coo_el.setter
    def k_coo_el(self, _):
        """ Sets the carbon monoxide - oxygen  gas-phase reaction parameter.
        """
        self.__k_coo_el = "k.COO.el"

    @k_coo_el.deleter
    def k_coo_el(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_o_ads(self):
        """ Gets the oxygen adsorption order.
        """
        return cp.deepcopy(self.__o_o_ads)

    @o_o_ads.setter
    def o_o_ads(self, _):
        """ Sets the oxygen adsorption parameter.
        """
        self.__o_o_ads = 2

    @o_o_ads.deleter
    def o_o_ads(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_o_des(self):
        """ Gets the oxygen desorption order.
        """
        return cp.deepcopy(self.__o_o_des)

    @o_o_des.setter
    def o_o_des(self, _):
        """ Sets the oxygen desorption parameter.
        """
        self.__o_o_des = 2

    @o_o_des.deleter
    def o_o_des(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_o_dif(self):
        """ Gets the oxygen diffusion order.
        """
        return cp.deepcopy(self.__o_o_dif)

    @o_o_dif.setter
    def o_o_dif(self, _):
        """ Sets the oxygen diffusion parameter.
        """
        self.__o_o_dif = 2

    @o_o_dif.deleter
    def o_o_dif(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_co_ads(self):
        """ Gets the carbon monoxide adsorption order.
        """
        return cp.deepcopy(self.__o_co_ads)

    @o_co_ads.setter
    def o_co_ads(self, _):
        """ Sets the carbon monoxide adsorption parameter.
        """
        self.__o_co_ads = 1

    @o_co_ads.deleter
    def o_co_ads(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_co_des(self):
        """ Gets the carbon monoxide desorption order.
        """
        return cp.deepcopy(self.__o_co_des)

    @o_co_des.setter
    def o_co_des(self, _):
        """ Sets the carbon monoxide desorption parameter.
        """
        self.__o_co_des = 1

    @o_co_des.deleter
    def o_co_des(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_co_dif(self):
        """ Gets the carbon monoxide diffusion order.
        """
        return cp.deepcopy(self.__o_co_dif)

    @o_co_dif.setter
    def o_co_dif(self, _):
        """ Sets the carbon monoxide diffusion parameter.
        """
        self.__o_co_dif = 2

    @o_co_dif.deleter
    def o_co_dif(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_coo_lh(self):
        """ Gets the carbon monoxide - oxygen reaction order.
        """
        return cp.deepcopy(self.__o_coo_lh)

    @o_coo_lh.setter
    def o_coo_lh(self, _):
        """ Sets the carbon monoxide - oxygen reaction parameter.
        """
        self.__o_coo_lh = 1

    @o_coo_lh.deleter
    def o_coo_lh(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_coo_el(self):
        """ Gets the carbon monoxide - oxygen  gas-phase reaction order.
        """
        return cp.deepcopy(self.__o_coo_el)

    @o_coo_el.setter
    def o_coo_el(self, _):
        """ Sets the carbon monoxide - oxygen  gas-phase reaction parameter.
        """
        self.__o_coo_el = 1

    @o_coo_el.deleter
    def o_coo_el(self):
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

        # Verify that the states allow a string representation and are unique
        if not isinstance(states, Iterable):
            raise ValueError(f"The states variable must be an iterable object, current type: {type(states)}")

        # Convert the iterable to a tuple of strings.
        tmp_states = tuple(map(str, states))

        # Check that the elements are unique.
        if not len(tmp_states) == len(tuple(set(states))):
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
                if len(state1[0]) >= process_orders[key1] and len(state1[1][key1]) > 0:
                    state0[1][key1].extend([cp.deepcopy(state1[0]) for _ in state1[1][key1]])

        def get_involved_states():
            """ Gets the states that are involved in the calculation.

                :return involved_states0: A list with the states that are
                involved in the calculation.
            """

            # Auxiliary variables.
            involved_orders0 = []
            involved_states0 = []

            # ------------------------------------------------------------------
            # Get the orders involved.
            # ------------------------------------------------------------------

            # Get the orders of the processes.
            orders0 = self._get_orders()

            # Get the order of the states for which the gains/decays will happen.
            length_states0 = tuple(set(tuple(map(len, lowest_states))))

            # Get the orders of the states involved in the calculations.
            for order0 in orders0:
                # Get the specific order of each state.
                for length_state0 in length_states0:
                    # Use the criteria to add the states.
                    involved_orders0.append(length_state0 + order0 - 1)

            # Filter the states, so that they are all unique.
            involved_orders0 = [order1 if order1 < self.number_of_sites else self.number_of_sites for order1 in involved_orders0]
            involved_orders0 = tuple(involved_orders0)
            involved_orders0 = set(involved_orders0)
            involved_orders0 = tuple(involved_orders0)

            # Get ALL the involved states.
            for order0_1 in involved_orders0:
                involved_states0.extend(self._get_states(order0_1))

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
            for i in range(2, order + 1):
                for involved_state0 in self._get_states(i):
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
                for i, process0 in enumerate(process_functions):
                    decay_dictionary_state0[keys[i]] = process0[2](state0)

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

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Get the zeroth order equations, these will serve as the basis.
        if order == 0 or order < self.number_of_sites:
            self.get_0th_order_equations(print_equations)

            # Get the equations to the lowest order.
            if order == 0 or order == 1:
                return

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

        # TODO: Work in the order reducer.

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

            # Add the equation to the list.
            # self.equations.append(cp.deepcopy(low_state))

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

    def _get_contracted_state(self, states, entry=-1):
        """ From a list of states, it returns the completely contracted state.
            For this to happen, the list of states must contain as much states
            as there arr possible number of states and the indexes of ALL the
            states must be the same.

            :param states: The list of states that are to be contracted.

            :param entry: The entry of the list to be contracted. It must be a
            number between zero and the length of one of the states, or a
            negative number, i.e., an index that indicates the index of the
            array to contracted.

            :return contracted_state: The contracted state. If it is empty, it
            means that the state cannot be contracted.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def validate_states(entry0):
            """ Checks that the states are valid to perform an index
                contraction operation.

                :param entry0: The entry0 of the list to be contracted. It must
                be a number between zero and the length of one of the states, or
                a negative number, i.e., an index that indicates the index of
                the array to contracted.
            """

            # Verify it is a list of valid states.
            for state0 in states:
                self._validate_state(state0)

            # Verify that the list contains as much states as possible site states.
            if not len(states) == len(self.states):
                raise ValueError("When a state is to be contracted, there must be "
                                 "as many states as particles. Number of requested = "
                                 f"{len(states)}, possible site states = {len(self.states)}")

            # Verify that the index is within the limits.
            while entry0 < 0:
                entry0 += len(states[0])

            # Verify that the index is
            if entry 0 > len(state0[0])



        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Validate the states
        validate_states()

        # Auxiliary variables.
        states_indexes = []
        states_particles = []

        # ----------------------------------------------------------------------
        # Validate the indexes for contraction.
        # ----------------------------------------------------------------------

    def _get_keys(self):
        """ Returns a tuple that contains the strings that represent the
            processes that serve as keys to the dictionaries of processes.

            :return keys: The strings that represent the processes that serve as
            keys to the dictionaries of processes.
        """

        # Get the process functions.
        functions, _ = self._get_process_functions()

        # The keys are in the second column.
        keys = tuple([key[1] for key in functions])

        return keys

    def _get_numbering(self, state):
        """ Returns a 2-tuple that contains the possible numbering schemes that
            the given state can take. The numbering scheme used is that of
            incremental consecutive numbering.

            :param state: The state to be numbered.

            :return: The list of possible numbered states in the given scheme.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def validate_unnumbered_state(state0):
            """ Validates that the length of the state is consistent.

                :param state0: The state to be validated.
            """

            # Validate the state.
            if 0 == len(state0) or len(state0) > self.number_of_sites:
                raise ValueError(f"The state must contain at least one site and at most {self.number_of_sites}. "
                                 f"Requested state sites: {len(state0)}")

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Validate the state.
        validate_unnumbered_state(state)

        # Auxiliary variables.
        all_states = []

        # Explore all the possibilities.
        for i in range(self.number_of_sites):
            # Only attempt if there are enough sites.
            if i + len(state) > self.number_of_sites:
                break

            # Make a deep copy of the state.
            tmp_state = cp.deepcopy(state)

            # Get the numbering list.
            tmp_list = list(range(i, i+len(state)))

            # Add the state to the list of possible states.
            all_states.append(tuple([(tmp_state[j], x + 1) for j, x in enumerate(tmp_list)]))

        return all_states

    def _get_orders(self):
        """ Gets the order of the states that are involved in the equations.
        """

        # Get the orders of the processes.
        tmp_list = tuple(x[0] for x in self._get_process_functions()[0])

        # Make the values unique.
        tmp_list = tuple(set(tmp_list))

        # Validate that the values are positive and less than the lattice site number.
        if any(map(lambda x:  0 >= x or x > self.number_of_sites, tmp_list)):
            raise ValueError("The orders in the equations must be greater than 0 "
                             f" and less than or equal to {self.number_of_sites}. ")

        return tmp_list

    def _get_process_functions(self):
        """ Gets the a 3-tuple with the process functions associated with the
            system, their minimum order and the string representation of the
            rate constant.

            :return process_functions: A 2-tuple that contains an N-tuple of
            3-tuples with the process functions associated with the system,
            their minimum order and the string representation of the rate
            constant; the second entry is an dictionary of the rate constant
            with the minimum order of the state for the process to happen.

            The tuples are in the order:

                process_functions: (process_order, rate_string, function_pointer)
                process_dictionary = (rate_string, process_order)
        """

        # Create the process function tuple.
        process_functions = (
            (self.o_o_ads, self.k_o_ads, self._oxygen_adsorb),
            (self.o_o_des, self.k_o_des, self._oxygen_desorb),
            (self.o_o_dif, self.k_o_dif, self._oxygen_diffusion),
            (self.o_co_ads, self.k_co_ads, self._carbon_monoxide_adsorb),
            (self.o_co_des, self.k_co_des, self._carbon_monoxide_desorb),
            (self.o_co_dif, self.k_co_dif, self._carbon_monoxide_diffusion),
            (self.o_coo_lh, self.k_coo_lh, self._lh_reaction),
            (self.o_coo_er, self.k_coo_er, self._er_reaction)
        )

        # Create the process dictionary.
        process_dictionary = ((process[1], process[0]) for process in process_functions)
        process_dictionary = dict(process_dictionary)

        return process_functions, process_dictionary

    def _get_states(self, order=1):
        """ Given the order, it returns a tuple of ALL the possible combinations
            of the system variables, i.e., all the possible combinations of the
            variables in N slots, where N=order; non-numbered.

            :param order: The order of the requested possible variables.

            :return all_states: All the variables of the given order.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def validate_order():
            """ Validates that the order to which the equations are requested
                are valid.
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

        # Obtain all the possible states.
        all_states = [cp.deepcopy(self.states) for _ in range(order)]
        all_states = product(*all_states)

        # ----------------------------------------------------------------------
        # Number the states properly.
        # ----------------------------------------------------------------------

        # Get a copy of ALL the possible states.
        tmp_states = cp.deepcopy(all_states)

        # Number the states.
        all_states = []
        for state in tmp_states:
            all_states.extend(self._get_numbering(state))

        # Turn the states into a tuple.
        all_states = tuple(all_states)

        return all_states

    def _get_state1_in_state2(self, state1, state2):
        """ Determines if state 1 is a sub-state of state 2. Order matters in
            this case.

            :param state1: The state that is evaluated to be a sub-state.

            :param state2: The state where state1 is going to be searched.

            :return: True if state1 is a substate, or the same, as state2.
            False, otherwise.
        """

        # Validate the states.
        self._validate_state(state1)
        self._validate_state(state2)

        # Determine if they are the same.
        is_sub_state = state1 == state2

        # If the states are equal no need to continue.
        if is_sub_state:
            return is_sub_state

        # If the state1 is longer than state 2, the state cannot be a
        # subtate of state2.
        if len(state1) > len(state2):
            return False

        # Otherwise, ALL the entries in state1 must be in state2.
        for i, entry in enumerate(state1):
            # The first entry can determine the state.
            if i == 0:
                is_sub_state = any(map(lambda x: entry == x, state2))
                continue

            # If an entry is not in the state, it cannot be a substate.
            if not is_sub_state:
                break

            # Check all the entries.
            is_sub_state = is_sub_state and any(map(lambda x: entry == x, state2))

        return is_sub_state

    # --------------------------------------------------------------------------
    # Carbon monoxide exclusive methods.
    # --------------------------------------------------------------------------

    def _carbon_monoxide_adsorb(self, initial_state):
        """
            Given a state, returns the list of states that are generated by the
            adsorption of carbon monoxide. Adsorption of carbon monoxide
            requires just a single empty site for the process to happen.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of carbon monoxide adsorption.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def append_to_final_states(state0):
            """ Appends the given state to the final state list, with the
                correct format.

                :param state0: The state to be appended.
            """
            # Get the final state.
            final_state = tuple((state1, tmp_numbering[j]) for j, state1 in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the state to the final state.
            final_states.append(final_state)

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Auxiliary variables.
        final_states = []

        # Check that the length of the state is consistent.
        if len(initial_state) < self.o_co_ads:
            return final_states

        # ----------------------------------------------------------------------
        # Get the list of states and the numbering.
        # ----------------------------------------------------------------------

        # Turn the 2-tuples in two lists.
        tmp_array = np.array(list(map(list, initial_state)))

        # List of unnumbered states.
        tmp_state = tmp_array[:, 0]

        # List of the numbering.
        tmp_numbering = np.array([int(x) for x in tmp_array[:, 1]])

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Check all the possible sites.
        for i, state in enumerate(tmp_state):
            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(tmp_state)

            # Adsorption is only possible if the site is empty.
            if state == 'E':
                tmp_state_0[i] = "CO"
                append_to_final_states(tmp_state_0)

        return final_states

    def _carbon_monoxide_desorb(self, initial_state):
        """ Given a state, returns the list of states that are generated by the
            desorption of carbon monoxide. Desorption of carbon monoxide
            requires just a single site with carbon monoxide for the process to
            happen.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of carbon monoxide desorption.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def append_to_final_states(state0):
            """ Appends the given state to the final state list, with the
                correct format.

                :param state0: The state to be appended.
            """

            # Get the final state.
            final_state = tuple((state1, tmp_numbering[j]) for j, state1 in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the state to the final state.
            final_states.append(final_state)

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Auxiliary variables.
        final_states = []

        # Check that the length of the state is consistent.
        if len(initial_state) < self.o_co_des:
            return final_states

        # ----------------------------------------------------------------------
        # Get the list of states and the numbering.
        # ----------------------------------------------------------------------

        # Turn the 2-tuples in two lists.
        tmp_array = np.array(list(map(list, initial_state)))

        # List of unnumbered states.
        tmp_state = tmp_array[:, 0]

        # List of the numbering.
        tmp_numbering = np.array([int(x) for x in tmp_array[:, 1]])

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Check all the possible sites.
        for i, state in enumerate(tmp_state):
            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(tmp_state)

            # Desorption is only possible if a carbon monoxide is present.
            if state == 'CO':
                tmp_state_0[i] = "E"
                append_to_final_states(tmp_state_0)

        return final_states

    def _carbon_monoxide_diffusion(self, initial_state):
        """
            Given a state, returns the list of states that are generated by the
            diffusion of carbon monoxide to adjacent sites. Adsorption of carbon
            monoxide requires two sites, one empty site and an adjacent site
            with carbon monoxide.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of carbon monoxide diffusion.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def append_to_final_states(state0):
            """ Appends the given state to the final state list, with the
                correct format.

                :param state0: The state to be appended.
            """

            # Get the final state.
            final_state = tuple((state1, tmp_numbering[j]) for j, state1 in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the state to the final state.
            final_states.append(final_state)

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Auxiliary variables.
        final_states = []

        # Check that the length of the state is consistent.
        if len(initial_state) < self.o_co_dif:
            return final_states

        # ----------------------------------------------------------------------
        # Get the list of states and the numbering.
        # ----------------------------------------------------------------------

        # Turn the 2-tuples in two lists.
        tmp_array = np.array(list(map(list, initial_state)))

        # List of unnumbered states.
        tmp_state = tmp_array[:, 0]

        # List of the numbering.
        tmp_numbering = np.array([int(x) for x in tmp_array[:, 1]])

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Check all the possible sites.
        for i, state in enumerate(tmp_state):
            # No need to traverse the whole state.
            if i == len(tmp_state) - 1:
                break

            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(tmp_state)

            # Diffusion is only possible if carbon monoxide is present and an adjacent site empty.
            if (tmp_state_0[i] == 'E' and tmp_state_0[i + 1] == 'CO') or (tmp_state_0[i] == 'CO' and tmp_state_0[i + 1] == 'E'):
                tmp_state_0[i], tmp_state_0[i + 1] = tmp_state_0[i + 1], tmp_state_0[i]
                append_to_final_states(tmp_state_0)

        return final_states

    # --------------------------------------------------------------------------
    # Oxygen exclusive methods.
    # --------------------------------------------------------------------------

    def _oxygen_adsorb(self, initial_state):
        """ Given a state, returns the list of states that are generated by the
            adsorption of oxygen. Adsorption of oxygen requires two sites, both
            sites must be neighbors and be empty for the process to happen.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of oxygen adsorption.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def append_to_final_states(state0):
            """ Appends the given state to the final state list, with the
                correct format.

                :param state0: The state to be appended.
            """

            # Get the final state.
            final_state = tuple((state1, tmp_numbering[j]) for j, state1 in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the state to the final state.
            final_states.append(final_state)

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Auxiliary variables.
        final_states = []

        # Check that the length of the state is consistent.
        if len(initial_state) < self.o_o_ads:
            return final_states

        # ----------------------------------------------------------------------
        # Get the list of states and the numbering.
        # ----------------------------------------------------------------------

        # No need to check states that have a length less than 2.
        if len(initial_state) < 2:
            return final_states

        # Turn the 2-tuples in two lists.
        tmp_array = np.array(list(map(list, initial_state)))

        # List of unnumbered states.
        tmp_state = tmp_array[:, 0]

        # List of the numbering.
        tmp_numbering = np.array([int(x) for x in tmp_array[:, 1]])

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Check all the possible sites.
        for i, state in enumerate(tmp_state):
            # Don't exceed the list length.
            if i == len(tmp_state) - 1:
                break

            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(tmp_state)

            # Adsorption is only possible if two adjacent empty sites are present.
            if tmp_state[i] == 'E' and tmp_state[i + 1] == 'E':
                tmp_state_0[i], tmp_state_0[i + 1] = "O", "O"
                append_to_final_states(tmp_state_0)

        return final_states

    def _oxygen_desorb(self, initial_state):
        """ Given a state, returns the list of states that are generated by the
            desorption of oxygen. Desorption of oxygen requires two sites, both
            sites must be neighbors and contain oxygen for the process to
            happen.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of oxygen desorption.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def append_to_final_states(state0):
            """ Appends the given state to the final state list, with the
                correct format.

                :param state0: The state to be appended.
            """

            # Get the final state.
            final_state = tuple((state1, tmp_numbering[j]) for j, state1 in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the state to the final state.
            final_states.append(final_state)

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Auxiliary variables.
        final_states = []

        # Check that the length of the state is consistent.
        if len(initial_state) < self.o_o_des:
            return final_states

        # ----------------------------------------------------------------------
        # Get the list of states and the numbering.
        # ----------------------------------------------------------------------

        # No need to check states that have a length less than 2.
        if len(initial_state) < 2:
            return final_states

        # Turn the 2-tuples in two lists.
        tmp_array = np.array(list(map(list, initial_state)))

        # List of unnumbered states.
        tmp_state = tmp_array[:, 0]

        # List of the numbering.
        tmp_numbering = np.array([int(x) for x in tmp_array[:, 1]])

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Check all the possible sites.
        for i, state in enumerate(tmp_state):
            # Don't exceed the list length.
            if i == len(tmp_state) - 1:
                break

            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(tmp_state)

            # Desorption is only possible if two adjacent oxygens are present.
            if tmp_state[i] == 'O' and tmp_state[i + 1] == 'O':
                tmp_state_0[i], tmp_state_0[i + 1] = "E", "E"
                append_to_final_states(tmp_state_0)

        return final_states

    def _oxygen_diffusion(self, initial_state):
        """ Given a state, returns the list of states that are generated by the
            diffusion of oxygen to adjacent sites. Adsorption of oxygen requires
            two sites, one empty site and an adjacent site with oxygen.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of oxygen diffusion.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def append_to_final_states(state0):
            """ Appends the given state to the final state list, with the
                correct format.

                :param state0: The state to be appended.
            """

            # Get the final state.
            final_state = tuple((state1, tmp_numbering[j]) for j, state1 in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the state to the final state.
            final_states.append(final_state)

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Auxiliary variables.
        final_states = []

        # Check that the length of the state is consistent.
        if len(initial_state) < self.o_o_dif:
            return final_states

        # ----------------------------------------------------------------------
        # Get the list of states and the numbering.
        # ----------------------------------------------------------------------

        # No need to check states that have a length less than 2.
        if len(initial_state) < 2:
            return final_states

        # Turn the 2-tuples in two lists.
        tmp_array = np.array(list(map(list, initial_state)))

        # List of unnumbered states.
        tmp_state = tmp_array[:, 0]

        # List of the numbering.
        tmp_numbering = np.array([int(x) for x in tmp_array[:, 1]])

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Check all the possible sites.
        for i, state in enumerate(tmp_state):
            # No need to traverse the whole state.
            if i == len(tmp_state) - 1:
                break

            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(tmp_state)

            # Oxygen diffusion is only possible if oxygen is present and an adjacent site empty.
            if (tmp_state_0[i] == 'E' and tmp_state_0[i + 1] == 'O') or (tmp_state_0[i] == 'O' and tmp_state_0[i + 1] == 'E'):
                tmp_state_0[i], tmp_state_0[i + 1] = tmp_state_0[i + 1], tmp_state_0[i]
                append_to_final_states(tmp_state_0)

        return final_states

    # --------------------------------------------------------------------------
    # Reaction Carbon monoxide - oxygen, exclusive methods.
    # --------------------------------------------------------------------------

    def _lh_reaction(self, initial_state):
        """ Given a state, returns the list of states that are generated by a
            Langmuir-Hinshelwood type of reaction. This reaction requires two
            adjacent sites, one with carbon monoxide and the other one with
            oxygen.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of a Langmuir-Hinashelwood reaction.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def append_to_final_states(state0):
            """ Appends the given state to the final state list, with the
                correct format.

                :param state0: The state to be appended.
            """

            # Get the final state.
            final_state = tuple((state1, tmp_numbering[j]) for j, state1 in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the state to the final state.
            final_states.append(final_state)

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Auxiliary variables.
        final_states = []

        # Check that the length of the state is consistent.
        if len(initial_state) < self.o_coo_lh:
            return final_states

        # ----------------------------------------------------------------------
        # Get the list of states and the numbering.
        # ----------------------------------------------------------------------

        # No need to check states that have a length less than 2.
        if len(initial_state) < 2:
            return final_states

        # Turn the 2-tuples in two lists.
        tmp_array = np.array(list(map(list, initial_state)))

        # List of unnumbered states.
        tmp_state = tmp_array[:, 0]

        # List of the numbering.
        tmp_numbering = np.array([int(x) for x in tmp_array[:, 1]])

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Check all the possible sites.
        for i, state in enumerate(tmp_state):
            # No need to traverse the whole state.
            if i == len(tmp_state) - 1:
                break

            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(tmp_state)

            # A Langmuir-Hinshelwood reaction is only possible if there are two
            # adjacent carbon monoxide and oxygen atoms.
            if (tmp_state_0[i] == 'CO' and tmp_state_0[i + 1] == 'O') or (tmp_state_0[i] == 'O' and tmp_state_0[i + 1] == 'CO'):
                tmp_state_0[i], tmp_state_0[i + 1] = "E", "E"
                append_to_final_states(tmp_state_0)

        return final_states

    def _er_reaction(self, initial_state):
        """ Given a state, returns the list of states that are generated by a an
            Elay-Rideal type of reaction. This reaction requires a single site
            with oxygen.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of an Elay-Rideal reaction.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def append_to_final_states(state0):
            """ Appends the given state to the final state list, with the
                correct format.

                :param state0: The state to be appended.
            """

            # Get the final state.
            final_state = tuple((state1, tmp_numbering[j]) for j, state1 in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the state to the final state.
            final_states.append(final_state)

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Auxiliary variables.
        final_states = []

        # Check that the length of the state is consistent.
        if len(initial_state) < self.o_coo_er:
            return final_states

        # ----------------------------------------------------------------------
        # Get the list of states and the numbering.
        # ----------------------------------------------------------------------

        # Turn the 2-tuples in two lists.
        tmp_array = np.array(list(map(list, initial_state)))

        # List of unnumbered states.
        tmp_state = tmp_array[:, 0]

        # List of the numbering.
        tmp_numbering = np.array([int(x) for x in tmp_array[:, 1]])

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Check all the possible sites.
        for i, state in enumerate(tmp_state):
            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(tmp_state)

            # An Elay-Rideal reaction is only possible if oxygen present.
            if tmp_state_0[i] == 'O':
                tmp_state_0[i] = "E"
                append_to_final_states(tmp_state_0)

        return final_states

    # --------------------------------------------------------------------------
    # Validation methods.
    # --------------------------------------------------------------------------

    def _validate_state(self, state):
        """ Validates that the state is a collection of 2-tuples and that
            the tuples have the form ("particle", "site"); where "particle" is in
            the set self.states and "site" is an integer in the range
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
                            f"a tuple. Types of state: {tmp_types}")

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

    def __init__(self):
        """ Creates a MathematicaGenerator object and initializes its
            properties.
        """

        # ----------------------------------------------------------------------
        # Define order of the process; i.e., minimum length of the state for the
        # process to happen.
        # ----------------------------------------------------------------------

        # Represents the minimum length of the state for oxygen related
        # processes.
        self.o_o_ads = 2
        self.o_o_des = 2
        self.o_o_dif = 2

        # Represents the minimum length of the state for carbon monoxide related
        # processes.
        self.o_co_ads = 1
        self.o_co_des = 1
        self.o_co_dif = 2

        # Represents the minimum length of the state for carbon monoxide -
        # oxygen reaction related processes.
        self.o_coo_lh = 2
        self.o_coo_er = 1

        # ----------------------------------------------------------------------
        # Define the default model parameters.
        # ----------------------------------------------------------------------

        # Define the maximum number of sites.
        self.number_of_sites = 3

        # Define the possible unique states each site of the system can take.
        self.states = ["CO", "O", "E"]

        # Array where the equations are saved.
        self.equations = []

        # ----------------------------------------------------------------------
        # Define the strings for the constant names.
        # ----------------------------------------------------------------------

        # String that represents the rate constants of oxygen related processes.
        self.k_o_ads = "k.O.ads"
        self.k_o_des = "k.O.des",
        self.k_o_dif = "k.O.Dif"

        # String that represents the rate constants of carbon monoxide related
        # processes.
        self.k_co_ads = "k.CO.ads"
        self.k_co_des = "k.CO.des"
        self.k_co_dif = "k.CO.dif"

        # String that represents the rate constants of carbon monoxide - oxygen
        # reaction related processes.
        self.k_coo_lh = "k.COO.lh"
        self.k_coo_er = "k.COO.er"

    # --------------------------------------------------------------------------
    # Dunder Methods.
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # Dunder Variables.
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


if __name__ == "__main__":
    # Test state:
    state_list0 = [(('E', 1), ('E', 2), ('E', 3)),
                   (('CO', 1), ('E', 2), ('E', 3)),
                   (('O', 1), ('E', 2), ('E', 3))
                   ]


    # Create the equation generator.
    tmp = EquationGenerator()

    # for i, state213 in enumerate(list_of_states):
    #     print(i, state213)

    tmp._get_contracted_state(state_list0, 2)


    # tmp.get_nth_order_equations(order=2, print_equations=False)
    #
    # tmp.generate_equations(gather_by_state=True, format_string="mathematica", order=1, save_file_name="tmp.txt")
