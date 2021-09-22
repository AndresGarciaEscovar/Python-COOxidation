""" Writes the equations for carbon monoxide - oxygen associative reaction on
    ruthenium (111).
"""

# Imports: General.
import copy as cp
import numpy as np

# Imports: Inherited and auxiliary user-defined classes.
from .mathematica_generator import EquationGenerator


class COOxidationEquationGenerator(EquationGenerator):
    """ Writes the equations in various formats, to different orders for the
        carbon monoxide - oxygen associative reaction on ruthenium (111);
        J. Chem. Phys. 143, 204702 (2015). https://doi.org/10.1063/1.4936354

        Class parameters:

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

        :param self.o_co_dif: An integer that represents the minimum length of
        the state for carbon monoxide diffusion.

        :param self.o_o_ads: An integer that represents the minimum length of
        the state for oxygen adsoprtion.

        :param self.o_o_des: An integer tnat represents the minimum length of
        the state for oxygen desoprtion.

        :param self.o_o_dif: An integer that represents the minimum length of
        the state for oxygen diffusion.

        Inherited parameters:

        :param self.sites: The maximum number of sites the system has.

        :param  self.states: The UNIQUE states in which each side can be in.
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
    def k_o_ads(self, _):
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
    def k_coo_er(self):
        """ Gets the carbon monoxide - oxygen  gas-phase reaction order.
        """
        return cp.deepcopy(self.__k_coo_er)

    @k_coo_er.setter
    def k_coo_er(self, _):
        """ Sets the carbon monoxide - oxygen  gas-phase reaction parameter.
        """
        self.__k_coo_er = "k.COO.el"

    @k_coo_er.deleter
    def k_coo_er(self):
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
    def o_coo_er(self):
        """ Gets the carbon monoxide - oxygen  gas-phase reaction order.
        """
        return cp.deepcopy(self.__o_coo_er)

    @o_coo_er.setter
    def o_coo_er(self, _):
        """ Sets the carbon monoxide - oxygen  gas-phase reaction parameter.
        """
        self.__o_coo_er = 1

    @o_coo_er.deleter
    def o_coo_er(self):
        pass

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Public Interface.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Private Interface.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Get Methods.
    # --------------------------------------------------------------------------

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
            (self.o_o_ads, self.k_o_ads, self._oxygen_adsorb,),
            (self.o_o_des, self.k_o_des, self._oxygen_desorb,),
            (self.o_o_dif, self.k_o_dif, self._oxygen_diffusion,),
            (self.o_co_ads, self.k_co_ads, self._carbon_monoxide_adsorb,),
            (self.o_co_des, self.k_co_des, self._carbon_monoxide_desorb,),
            (self.o_co_dif, self.k_co_dif, self._carbon_monoxide_diffusion,),
            (self.o_coo_lh, self.k_coo_lh, self._lh_reaction,),
            (self.o_coo_er, self.k_coo_er, self._er_reaction,)
        )

        return process_information

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

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def validate_unnumbered_state(state0):
            """ Validates that the length of the state is consistent.

                :param state0: The state to be validated.
            """

            # Validate the state.
            if len(state0) == 0 or len(state0) > self.number_of_sites:
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
            tmp_list = list(range(i, i + len(state)))

            # Add the state to the list of possible states.
            all_states.append(tuple((tmp_state[j], x + 1) for j, x in enumerate(tmp_list)))

        return all_states

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

    # --------------------------------------------------------------------------
    # Carbon monoxide exclusive methods.
    # --------------------------------------------------------------------------

    def _carbon_monoxide_adsorb(self, initial_state):
        """ Given a state, returns the list of states that are generated by the
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
            final_state = tuple((cp.deepcopy(entry), cp.deepcopy(tmp_indexes[j])) for j, entry in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the list of possible final states.
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

        # Get the numpy decomposition of the state.
        tmp_particles, tmp_indexes = self._get_state_elements(initial_state)

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Check all the possible sites.
        for i, state in enumerate(tmp_particles):
            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(list(tmp_particles))

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
            final_state = tuple((cp.deepcopy(entry), cp.deepcopy(tmp_indexes[j])) for j, entry in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the list of possible final states.
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

        # Get the numpy decomposition of the state.
        tmp_particles, tmp_indexes = self._get_state_elements(initial_state)

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Check all the possible sites.
        for i, state in enumerate(tmp_particles):
            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(list(tmp_particles))

            # Desorption is only possible if the site contains a carbon monoxide.
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
            final_state = tuple((cp.deepcopy(entry), cp.deepcopy(tmp_indexes[j])) for j, entry in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the list of possible final states.
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

        # Get the numpy decomposition of the state.
        tmp_particles, tmp_indexes = self._get_state_elements(initial_state)

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Check all the possible sites.
        for i, state in enumerate(list(tmp_particles)):
            # Get out of the loop before the last element.
            if i == len(tmp_particles) - 1:
                break

            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(list(tmp_particles))

            # Can only diffuse to adjacent sites.
            indexes_valid = tmp_indexes[i] + 1 == tmp_indexes[i + 1]

            # Can only happen with an empty site and a site with carbon monoxide.
            particles_valid = tmp_particles[i] == "E" and tmp_particles[i + 1] == "CO"
            particles_valid = particles_valid or (tmp_particles[i] == "CO" and tmp_particles[i + 1] == "E")

            # Check to see if carbon monoxide diffusion happens.
            if particles_valid and indexes_valid:
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
            final_state = tuple((cp.deepcopy(entry), cp.deepcopy(tmp_indexes[j])) for j, entry in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the list of possible final states.
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

        # Get the numpy decomposition of the state.
        tmp_particles, tmp_indexes = self._get_state_elements(initial_state)

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Check all the possible sites.
        for i, _ in enumerate(tmp_particles):
            # Get out of the loop before the last element.
            if i == len(tmp_particles) - 1:
                break

            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(list(tmp_particles))

            # Can only associate with adjacent sites.
            indexes_valid = tmp_indexes[i] + 1 == tmp_indexes[i + 1]

            # Can only happen with two adjacent oxygens.
            particles_valid = tmp_particles[i] == "E" and tmp_particles[i + 1] == "E"

            # Check if oxygen adsorption can happen.
            if particles_valid and indexes_valid:
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
            final_state = tuple((cp.deepcopy(entry), cp.deepcopy(tmp_indexes[j])) for j, entry in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the list of possible final states.
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

        # Get the numpy decomposition of the state.
        tmp_particles, tmp_indexes = self._get_state_elements(initial_state)

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Check all the possible sites.
        for i, _ in enumerate(tmp_particles):
            # Get out of the loop before the last element.
            if i == len(tmp_particles) - 1:
                break

            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(list(tmp_particles))

            # Can only associate with adjacent sites.
            indexes_valid = tmp_indexes[i] + 1 == tmp_indexes[i + 1]

            # Can only happen with two adjacent oxygens.
            particles_valid = tmp_particles[i] == "O" and tmp_particles[i + 1] == "O"

            # Check if oxygen desorption can happen.
            if particles_valid and indexes_valid:
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
            final_state = tuple((cp.deepcopy(entry), cp.deepcopy(tmp_indexes[j])) for j, entry in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the list of possible final states.
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

        # Get the numpy decomposition of the state.
        tmp_particles, tmp_indexes = self._get_state_elements(initial_state)

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Check all the possible sites.
        for i, _ in enumerate(tmp_particles):
            # Get out of the loop before the last element.
            if i == len(tmp_particles) - 1:
                break

            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(list(tmp_particles))

            # Can only diffuse to adjacent sites.
            indexes_valid = tmp_indexes[i] + 1 == tmp_indexes[i + 1]

            # Can only happen with an empty site and a site with oxygen.
            particles_valid = tmp_particles[i] == "E" and tmp_particles[i + 1] == "O"
            particles_valid = particles_valid or (tmp_particles[i] == "O" and tmp_particles[i + 1] == "E")

            # Check if oxygen diffusion can happen.
            if particles_valid and indexes_valid:
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
            to the process of a Langmuir-Hinshelwood reaction.
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
            final_state = tuple((cp.deepcopy(entry), cp.deepcopy(tmp_indexes[j])) for j, entry in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the list of possible final states.
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

        # Get the numpy decomposition of the state.
        tmp_particles, tmp_indexes = self._get_state_elements(initial_state)

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Check all the possible sites.
        for i, _ in enumerate(tmp_particles):
            # Get out of the loop before the last element.
            if i == len(tmp_particles) - 1:
                break

            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(list(tmp_particles))

            # Can only react with adjacent sites.
            indexes_valid = tmp_indexes[i] + 1 == tmp_indexes[i + 1]

            # Can only happen with adjacent carbon moxide and oxygen filled sites.
            particles_valid = tmp_particles[i] == "CO" and tmp_particles[i + 1] == "O"
            particles_valid = particles_valid or (tmp_particles[i] == "O" and tmp_particles[i + 1] == "CO")

            # Check if the Langmuir-Hinshelwood reaction can happen.
            if particles_valid and indexes_valid:
                tmp_state_0[i], tmp_state_0[i + 1] = "E", "E"
                append_to_final_states(tmp_state_0)

        return final_states

    def _er_reaction(self, initial_state):
        """ Given a state, returns the list of states that are generated by an
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
            final_state = tuple((cp.deepcopy(entry), cp.deepcopy(tmp_indexes[j])) for j, entry in enumerate(state0))

            # Validate the final state.
            self._validate_state(final_state)

            # Append to the list of possible final states.
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

        # Get the numpy decomposition of the state.
        tmp_particles, tmp_indexes = self._get_state_elements(initial_state)

        # ----------------------------------------------------------------------
        # Get the resulting states.
        # ----------------------------------------------------------------------

        # Check all the possible sites.
        for i, state in enumerate(tmp_particles):
            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(list(tmp_particles))

            # An Elay-Rideal reaction is only possible if the site contains oxygen.
            if state == 'O':
                tmp_state_0[i] = "E"
                append_to_final_states(tmp_state_0)

        return final_states

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Constructors and Dunder Methods.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Constructor.
    # --------------------------------------------------------------------------

    def __init__(self, sites=1, states=("E",)):
        """ Builds a class that writes the equations for the carbon monoxide -
            oxygen associative reaction on ruthenium (111);
            J. Chem. Phys. 143, 204702 (2015).

            Initializes the class with the standard parameters.
        """

        # Initialize the super class.
        super(COOxidationEquationGenerator, self).__init__(sites, ('CO', 'O', 'E'))

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