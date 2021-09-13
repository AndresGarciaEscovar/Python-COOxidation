""" Writes the equations for carbon monoxide - oxygen associative reaction on
    ruthenium (111); J. Chem. Phys. 143, 204702 (2015).
    https://doi.org/10.1063/1.4936354
"""

# Imports.
import copy as cp
from itertools import product
import numpy as np


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
    # Public Interface.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Get methods.
    # --------------------------------------------------------------------------

    def get_exact_equations(self):
        """ Gets the lowest order exact equations.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def get_resulting_state(state0):
            """ Gets a dictionary of resulting states from applying ALL the
                possible system operations to a given state.

                :param state0: The state to which the process will be applied.

                :return: A dictionary of the resulting states from applying ALL
                the possible system operations to the given state.
            """

            # Check that it is a valid state.
            self._validate_state(state0)

            # The resultant list of states when a process is applied.
            result_list = []

            # Attempt to apply all the processes to the given state.
            for j, process in enumerate(process_functions):
                # Apply the process to the given state.
                result_states = [] if not len(state0) == process[0] else process[2](state0)

                # Append the results to the list.
                result_list.append((process[1], result_states))

            # Create a dictionary
            result_list = tuple(result_list)

            return dict(result_list)

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Get the information of the process functions.
        process_functions = self._get_process_functions()

        # Get the lowest state variables, with the numbering.
        lowest_states = self._get_states(1)

        # Get ALL the states that are potentially involved in the calculation.
        involved_states = []
        for order in self._get_orders():
            involved_states.extend(self._get_states(order))

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # List where to store the results.
        resultant_states = []

        # Apply the functions to the ALL the states.
        for i, involved_state in enumerate(involved_states):
            # Store the results here.
            resultant_states.append([involved_state, get_resulting_state(involved_state)])

        # Get the dictionary keys from the first entry.
        reaction_constants = resultant_states[0][1].keys()

        # Loop through all the lowest order states.
        for i, low_state in enumerate(lowest_states):
            # Loop through all the possible outcomes of the processes for ALL the involved states.
            for j, resultant_state in enumerate(resultant_states):
                # If the low state is a substate of the state undergoing the process, it has the potential to
                # destroy the state.
                if self._get_state1_in_state2(low_state, resultant_state[0]):
                    # print(low_state, resultant_state[0])
                    continue

                # If the state is not a substate of the resultant state, it has the potential of creating the state.
                if self._get_state1_in_state2_indexes(low_state, resultant_state[0]):
                    print(low_state, resultant_state[0])

            print("")

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Private Interface.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Get methods.
    # --------------------------------------------------------------------------

    def _get_orders(self):
        """ Gets the order of the states that are involved in the equations.
        """

        # Oxygen related processes.
        tmp_list = [self.o_o_ads, self.o_o_des, self.o_o_dif]

        # Carbon monoxide processes.
        tmp_list.extend([self.o_co_ads, self.o_co_des, self.o_co_dif])

        # Carbon monoxide - oxygen reaction processes.
        tmp_list.extend([self.o_coo_lh, self.o_coo_er])

        # Gets the non-repeated elements of the list.
        tmp_list = tuple(set(tmp_list))

        # Validate that the values are positive and less than the given order.
        if any(map(lambda x:  0 >= x or x > self.sites, tmp_list)):
            raise ValueError("The orders in the equations must be greater than 0 "
                             f" and less than or equal to {self.sites}. ")

        return tmp_list

    def _get_process_functions(self):
        """ Gets the a 3-tuple with the process functions associated with the
            system, their minimum order and the string representation of the
            rate constant.

            :return process_functions: A list of 3-tuples with the process
            functions associated with the system, their minimum order and the
            string representation of the rate constant.

            The tuples are in the order:
            (process_order, rate_string, function_pointer)
        """
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

        # List of tuples with the associated rate of the process and the
        # pointer of the function of the specific process.

        return process_functions

    def _get_states(self, order=1):
        """ Given the order, it returns ALL the possible combinations of the
            system variables, i.e., all the possible combinations of the
            variables in n slots. Non-numbered.

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
            if order > self.sites:
                raise ValueError(f"The order parameter must less than or equal to {self.sites}.")

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

    def _get_state1_in_state2_indexes(self, state1, state2):
        """ Determines if indexes of the sites in state 1 are a subset of the
            indexes of the sites in state 2. Order matters in this case.

            :param state1: The state whose indexes are evaluated to be a
            sub-set.

            :param state2: The state where the indexes of state1 are going to be
            searched.

            :return: True if the indexes of state1 are a subset, or the same,
            as those of state2. False, otherwise.
        """

        # Validate the states.
        self._validate_state(state1)
        self._validate_state(state2)

        # Get the subindexes of state 1 and state 2.
        subindexes1 = tuple(x[1] for x in state1)
        subindexes2 = tuple(x[1] for x in state2)

        # Determine if they are the same.
        is_sub_state = subindexes1 == subindexes2

        # If the states are equal no need to continue.
        if is_sub_state:
            return is_sub_state

        # If the state1 is longer than state 2, the state cannot be a
        # subtate of state2.
        if len(subindexes1) > len(subindexes2):
            return False

        # Otherwise, ALL the entries in state1 must be in state2.
        for i, entry in enumerate(subindexes1):
            # The first entry can determine the state.
            if i == 0:
                is_sub_state = any(map(lambda x: entry == x, subindexes2))
                continue

            # If an entry is not in the state, it cannot be a substate.
            if not is_sub_state:
                break

            # Check all the entries.
            is_sub_state = is_sub_state and any(map(lambda x: entry == x, subindexes2))

        return is_sub_state

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
            if 0 == len(state0) or len(state0) > self.sites:
                raise ValueError(f"The state must contain at least one site and at most {self.sites}. "
                                 f"Requested state sites: {len(state0)}")

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Validate the state.
        validate_unnumbered_state(state)

        # Auxiliary variables.
        all_states = []

        # Explore all the possibilities.
        for i in range(self.sites):
            # Only attempt if there are enough sites.
            if i + len(state) > self.sites:
                break

            # Make a deep copy of the state.
            tmp_state = cp.deepcopy(state)

            # Get the numbering list.
            tmp_list = list(range(i, i+len(state)))

            # Add the state to the list of possible states.
            all_states.append([(tmp_state[j], x + 1) for j, x in enumerate(tmp_list)])

        return all_states

    def _get_is_site_in_state(self, site, state):
        """ Given a site and a state, it determines if the state is related to
            the site, i.e., the state has a numbering such that the number of
            the site appears.

            :param site: The site that is being searched.

            :param state: The state that is properly labeled.

            :return: True, if the site in the equation. False, otherwise.
        """

        # ----------------------------------------------------------------------
        # Auxiliary states.
        # ----------------------------------------------------------------------

        def validate_state():
            """ Validates if the state is valid, i.e., has the appropriate
                length.
            """

            # The length of the state must be greater than zero and, at most, self.sites.
            if len(state) == 0 or len(state) > self.sites:
                raise ValueError(f"The size of the state must be, at most, {self.sites}"
                                 f"Current size: {len(state)}")

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Validate the state.
        validate_state()

        # Get the numbering of the states.
        sites = [x[1] for x in state]

        return site in sites

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
            final_state =tuple((state1, tmp_numbering[j]) for j, state1 in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the state to the final state.
            final_states.append(final_state)

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Auxiliary variables.
        final_states = []

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

            # A Langmuir-Hinshelwood reaction is only possible if there are two adjacent carbon monoxide and oxygen atoms.
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
        if not 1 <= len(state) <= self.sites:
            raise ValueError(f"The current length of the state list is not valid, it must be"
                             f"in the range [1, {self.sites}].  Current legth: {len(state)}.")
        # ----------------------------------------------------------------------
        # Set the auxiliary variables.
        # ----------------------------------------------------------------------

        # Get the particles.
        state = np.array(state)

        # Get the particles and the numbering of the states.
        particles = list(state[:,0])
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
        if not 1 <= len(numbering) <= self.sites:
            raise ValueError(f"The current length of the numbering list is not valid, it must be"
                             f"in the range [1, {self.sites}].  Current length: {len(numbering)}.")

        # Check that the numbering of the state is greater than zero and less than or equal to the maximum number of
        # sites.
        if not all(map(lambda x: 0 < x <= self.sites, numbering)):
            raise ValueError(f"The numbering of the states must be in the range [1, {self.sites}]."
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
        self.sites = 3

        # Define the possible unique states each site of the system can take.
        self.states = ["CO", "O", "E"]

        # ----------------------------------------------------------------------
        # Define the strings for the constant names.
        # ----------------------------------------------------------------------

        # String that represents the rate constants of oxygen related processes.
        self.k_o_ads = "kOAds"
        self.k_o_des = "kODes"
        self.k_o_dif = "kODif"

        # String that represents the rate constants of carbon monoxide related
        # processes.
        self.k_co_ads = "kCOAds"
        self.k_co_des = "kCODes"
        self.k_co_dif = "kCODif"

        # String that represents the rate constants of carbon monoxide - oxygen
        # reaction related processes.
        self.k_coo_lh = "kCOOLH"
        self.k_coo_er = "kCOOER"

    # --------------------------------------------------------------------------
    # Dunder Methods.
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # Dunder Variables.
    # --------------------------------------------------------------------------

    __slots__ = ['k_o_ads', 'k_o_des', 'k_o_dif', 'k_co_ads', 'k_co_des', 'k_co_dif', 'k_coo_lh', 'k_coo_er',
                 'o_o_ads', 'o_o_des', 'o_o_dif', 'o_co_ads', 'o_co_des', 'o_co_dif', 'o_coo_lh', 'o_coo_er',
                 'sites', 'states']


if __name__ == "__main__":
    tmp = EquationGenerator()

    tmp.get_exact_equations()



