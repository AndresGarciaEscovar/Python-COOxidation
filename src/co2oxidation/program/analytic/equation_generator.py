""" Writes the equations for carbon monoxide - oxygen associative reaction on
    ruthenium (111).
"""

# ------------------------------------------------------------------------------
# Imports.
# ------------------------------------------------------------------------------

# Imports: General.
import copy as cp
import itertools
import os

# Imports: User-defined.
from coOxidation.Program.Analytic.Formatters.formatter_manager import FormatterManager

from coOxidation.Program.Analytic.Interfaces.generator import Generator

# ------------------------------------------------------------------------------
# Classes.
# ------------------------------------------------------------------------------


class EquationGenerator(Generator):
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

        :param self.constraints: The list where the constraint equations will be
         saved.

        :param self.equations: The list where the equations will be saved.

        :param self.sites: The maximum number of sites.

        :param  self.states: The UNIQUE states in which each site can be in.
    """

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Getters, Setters and Deleters.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    @property
    def k_o_ads(self) -> str:
        """ Gets the oxygen adsorption order.

            :return: The oxygen adsorption parameter.
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
    def k_o_des(self) -> str:
        """ Gets the oxygen desorption order.

            :return: The oxygen desorption parameter.
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
    def k_o_dif(self) -> str:
        """ Gets the oxygen diffusion order.

            :return: The oxygen diffusion parameter.
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
    def k_co_ads(self) -> str:
        """ Gets the carbon monoxide adsorption order.

            :return: The carbon monoxide adsorption parameter.
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
    def k_co_des(self) -> str:
        """ Gets the carbon monoxide desorption order.

            :return: The carbon monoxide desorption parameter.
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
    def k_co_dif(self) -> str:
        """ Gets the carbon monoxide diffusion order.

            :return: The carbon monoxide diffusion parameter.
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
    def k_coo_lh(self) -> str:
        """ Gets the carbon monoxide - oxygen reaction order.

            :return: The Langmuir-Hinshelwood reaction parameter.
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
    def k_coo_er(self) -> str:
        """ Gets the carbon monoxide - oxygen  gas-phase reaction order.

            :return: The Elay-Rideal reaction parameter.
        """
        return cp.deepcopy(self.__k_coo_er)

    @k_coo_er.setter
    def k_coo_er(self, _):
        """ Sets the carbon monoxide - oxygen  gas-phase reaction parameter.
        """
        self.__k_coo_er = "k.COO.er"

    @k_coo_er.deleter
    def k_coo_er(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_o_ads(self) -> int:
        """ Gets the oxygen adsorption order.

            :return: The oxygen adsorption order.
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
    def o_o_des(self) -> int:
        """ Gets the oxygen desorption order.

            :return: The oxygen desorption order.
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
    def o_o_dif(self) -> int:
        """ Gets the oxygen diffusion order.

            :return: The oxygen diffusion order.
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
    def o_co_ads(self) -> int:
        """ Gets the carbon monoxide adsorption order.

            :return: The oxygen adsorption order.
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
    def o_co_des(self) -> int:
        """ Gets the carbon monoxide desorption order.

            :return: The carbon monoxide desorption order.
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
    def o_co_dif(self) -> int:
        """ Gets the carbon monoxide diffusion order.

            :return: The carbon monoxide diffusion order.
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
    def o_coo_lh(self) -> int:
        """ Gets the carbon monoxide - oxygen reaction order.

            :return: The Langmuir-Hinshelwood reaction order.
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
    def o_coo_er(self) -> int:
        """ Gets the carbon monoxide - oxygen  gas-phase reaction order.

            :return: The Elay-Rideal reaction order.
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
    # Methods.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Get Methods.
    # --------------------------------------------------------------------------

    def get_associated_operations(self) -> tuple:
        """ Returns a tuple with the 3-tuples that contain the operations that
            can be applied to the states of the system, the order of the process
            and the string representation of the state.

            :return: A tuple with the 3-tuples that contain the operations that
             can be applied to the states of the system, the order of the
             process and the string representation of the state, in the order:
             (process order, process rate constant, pointer to function)
        """

        # Define the operations.
        process_information = (
            (self.o_o_ads, self.k_o_ads, self.process_oxygen_adsorb,),
            (self.o_o_des, self.k_o_des, self.process_oxygen_desorb,),
            (self.o_o_dif, self.k_o_dif, self.process_oxygen_diffusion,),
            (self.o_co_ads, self.k_co_ads, self.process_carbon_monoxide_adsorb,),
            (self.o_co_des, self.k_co_des, self.process_carbon_monoxide_desorb,),
            (self.o_co_dif, self.k_co_dif, self.process_carbon_monoxide_diffusion,),
            (self.o_coo_lh, self.k_coo_lh, self.process_reaction_lh_reaction,),
            (self.o_coo_er, self.k_coo_er, self.process_reaction_er_reaction,)
        )

        return process_information

    def get_constraints(self, states: list) -> None:
        """ Given a set of numbered states, it stores the constraints of the
            system in the self.constrainst variable, i.e., the probability
            identities, 1 = sum(x in X) P(x), P(x) = sum(y in Y) P(x,y),
            P(y) = sum(x in X) P(x,y), etc; with 0 <= P(x) <= 1 for x in X.


            :param states: ALL of the "left-hand" states of the system.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def expand_state(state0: tuple, site0: int) -> list:
            """ Returns the expanded state at the given site.

                :param state0: The state that is to be expanded.

                :param site0: The site at which it must be expanded.

                :return: The expanded state.
            """

            # List where the constraints will be placed.
            constraint_list0 = []

            # The index of the expanded states.
            index_to_add = state0[site0][1]
            index_to_add += 1 if site0 == -1 else -1

            # Normalize the index.
            site0 = 0 if site0 == 0 else len(state0)

            # Add the states with the indexes.
            for state0_0 in self.states:
                expanded_state0 = list(cp.deepcopy(state0))
                new_state = (state0_0, index_to_add,)
                expanded_state0.insert(site0, new_state)
                constraint_list0.append(tuple(expanded_state0))

            return constraint_list0

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        state_lengths = list(set(len(state) for state in states))
        state_lengths = sorted(state_lengths)[:-1]

        if len(state_lengths) == 0:
            return

        max_length = max(state_lengths)
        for state in states:
            if len(state) <= max_length:
                index = state[-1][1]
                constraint_list = expand_state(state, -1) if index < self.sites_number else expand_state(state, 0)
                self.constraints.append((state, constraint_list,))

    def get_numbering(self, state: tuple):
        """ Returns a tuple with the possible numbering a state of length N can
            have. The format of a SINGLE state must be given in the format:
            ( (particle_at_site1, numbering_scheme1),..., (particle_at_siteN,
            numbering_schemeN)).

            :param state: The state to be numbered.

            :return: The list of possible numbered states in the given format.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def validate_state(state0: tuple) -> None:
            """ From a list of

                :param state0: The state to be numbered.
            """

            if not (0 < len(state0) <= self.sites_number):
                raise ValueError(
                    "The state to be numbered must contain at least one element"
                    " and must not be longer the length of the lattice."
                )

            if not all(map(lambda x: x in self.states, state0)):
                raise TypeError("The types in the state are not all valid.")

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        validate_state(state)

        all_states = []
        for i in range(self.sites_number):
            if i + len(state) > self.sites_number:
                break

            tmp_state = cp.deepcopy(state)
            tmp_list = list(range(i, i + len(state)))
            all_states.append(tuple((tmp_state[j], x + 1) for j, x in enumerate(tmp_list)))

        return all_states

    # --------------------------------------------------------------------------
    # Process Carbon Monoxide Methods.
    # --------------------------------------------------------------------------

    def process_carbon_monoxide_adsorb(self, initial_state: tuple):
        """ Given a state, returns the list of states that are generated by the
            adsorption of carbon monoxide. Adsorption of carbon monoxide
            requires just a single empty site for the process to happen.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of carbon monoxide adsorption.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def append_to_final_states(state0: tuple):
            """ Appends the given state to the final state list, with the
                correct format.

                :param state0: The state to be appended.
            """

            final_state = tuple((cp.deepcopy(entry), cp.deepcopy(tmp_indexes[j])) for j, entry in enumerate(state0))
            final_states.append(final_state)

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        final_states = []
        if len(initial_state) < self.o_co_ads:
            return final_states

        # Get the list of states and the numbering.
        tmp_particles, tmp_indexes = self.get_state_elements(initial_state)
        for i, state in enumerate(tmp_particles):
            # Always make a copy of the initial state.
            tmp_state_0 = cp.deepcopy(list(tmp_particles))
            if state == 'E':
                tmp_state_0[i] = "CO"
                append_to_final_states(tuple(tmp_state_0))

        return final_states

    def process_carbon_monoxide_desorb(self, initial_state: tuple):
        """ Given a state, returns the list of states that are generated by the
            desorption of carbon monoxide. Desorption of carbon monoxide
            requires just a single site with carbon monoxide for the process to
            happen.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of carbon monoxide desorption.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def append_to_final_states(state0: tuple):
            """ Appends the given state to the final state list, with the
                correct format.

                :param state0: The state to be appended.
            """

            final_state = tuple((cp.deepcopy(entry), cp.deepcopy(tmp_indexes[j])) for j, entry in enumerate(state0))
            final_states.append(final_state)

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        final_states = []
        if len(initial_state) < self.o_co_des:
            return final_states

        # Get the list of states and the numbering.
        tmp_particles, tmp_indexes = self.get_state_elements(initial_state)

        # Get the resulting states.
        for i, state in enumerate(tmp_particles):
            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(list(tmp_particles))
            if state == 'CO':
                tmp_state_0[i] = "E"
                append_to_final_states(tuple(tmp_state_0))

        return final_states

    def process_carbon_monoxide_diffusion(self, initial_state: tuple):
        """
            Given a state, returns the list of states that are generated by the
            diffusion of carbon monoxide to adjacent sites. Adsorption of carbon
            monoxide requires two sites, one empty site and an adjacent site
            with carbon monoxide.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of carbon monoxide diffusion.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def append_to_final_states(state0: tuple):
            """ Appends the given state to the final state list, with the
                correct format.

                :param state0: The state to be appended.
            """

            final_state = tuple((cp.deepcopy(entry), cp.deepcopy(tmp_indexes[j])) for j, entry in enumerate(state0))
            final_states.append(final_state)

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        final_states = []
        if len(initial_state) < self.o_co_dif:
            return final_states

        # Get the list of states and the numbering.
        tmp_particles, tmp_indexes = self.get_state_elements(initial_state)

        # Get the resulting states.
        for i, state in enumerate(list(tmp_particles)):
            if i == len(tmp_particles) - 1:
                break

            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(list(tmp_particles))
            indexes_valid = tmp_indexes[i] + 1 == tmp_indexes[i + 1]
            particles_valid = tmp_particles[i] == "E" and tmp_particles[i + 1] == "CO"
            particles_valid = particles_valid or (tmp_particles[i] == "CO" and tmp_particles[i + 1] == "E")

            # Check to see if carbon monoxide diffusion happens.
            if particles_valid and indexes_valid:
                tmp_state_0[i], tmp_state_0[i + 1] = tmp_state_0[i + 1], tmp_state_0[i]
                append_to_final_states(tuple(tmp_state_0))

        return final_states

    # --------------------------------------------------------------------------
    # Process Oxygen Methods.
    # --------------------------------------------------------------------------

    def process_oxygen_adsorb(self, initial_state: tuple):
        """ Given a state, returns the list of states that are generated by the
            adsorption of oxygen. Adsorption of oxygen requires two sites, both
            sites must be neighbors and be empty for the process to happen.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of oxygen adsorption.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def append_to_final_states(state0: tuple):
            """ Appends the given state to the final state list, with the
                correct format.

                :param state0: The state to be appended.
            """

            final_state = tuple((cp.deepcopy(entry), cp.deepcopy(tmp_indexes[j])) for j, entry in enumerate(state0))
            final_states.append(final_state)

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        final_states = []
        if len(initial_state) < self.o_o_ads:
            return final_states

        # Get the list of states and the numbering.
        tmp_particles, tmp_indexes = self.get_state_elements(initial_state)
        for i, _ in enumerate(tmp_particles):
            if i == len(tmp_particles) - 1:
                break

            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(list(tmp_particles))
            indexes_valid = tmp_indexes[i] + 1 == tmp_indexes[i + 1]
            particles_valid = tmp_particles[i] == "E" and tmp_particles[i + 1] == "E"

            # Check if oxygen adsorption can happen.
            if particles_valid and indexes_valid:
                tmp_state_0[i], tmp_state_0[i + 1] = "O", "O"
                append_to_final_states(tuple(tmp_state_0))

        return final_states

    def process_oxygen_desorb(self, initial_state: tuple):
        """ Given a state, returns the list of states that are generated by the
            desorption of oxygen. Desorption of oxygen requires two sites, both
            sites must be neighbors and contain oxygen for the process to
            happen.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of oxygen desorption.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def append_to_final_states(state0: tuple):
            """ Appends the given state to the final state list, with the
                correct format.

                :param state0: The state to be appended.
            """

            final_state = tuple((cp.deepcopy(entry), cp.deepcopy(tmp_indexes[j])) for j, entry in enumerate(state0))
            final_states.append(final_state)

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        final_states = []
        if len(initial_state) < self.o_o_des:
            return final_states

        # Get the list of states and the numbering.
        tmp_particles, tmp_indexes = self.get_state_elements(initial_state)
        for i, _ in enumerate(tmp_particles):
            if i == len(tmp_particles) - 1:
                break

            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(list(tmp_particles))
            indexes_valid = tmp_indexes[i] + 1 == tmp_indexes[i + 1]
            particles_valid = tmp_particles[i] == "O" and tmp_particles[i + 1] == "O"
            if particles_valid and indexes_valid:
                tmp_state_0[i], tmp_state_0[i + 1] = "E", "E"
                append_to_final_states(tuple(tmp_state_0))

        return final_states

    def process_oxygen_diffusion(self, initial_state: tuple):
        """ Given a state, returns the list of states that are generated by the
            diffusion of oxygen to adjacent sites. Adsorption of oxygen requires
            two sites, one empty site and an adjacent site with oxygen.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of oxygen diffusion.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def append_to_final_states(state0: tuple):
            """ Appends the given state to the final state list, with the
                correct format.

                :param state0: The state to be appended.
            """

            final_state = tuple((cp.deepcopy(entry), cp.deepcopy(tmp_indexes[j])) for j, entry in enumerate(state0))
            final_states.append(final_state)

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        final_states = []
        if len(initial_state) < self.o_o_dif:
            return final_states

        # Get the list of states and the numbering.
        tmp_particles, tmp_indexes = self.get_state_elements(initial_state)
        for i, _ in enumerate(tmp_particles):
            if i == len(tmp_particles) - 1:
                break

            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(list(tmp_particles))
            indexes_valid = tmp_indexes[i] + 1 == tmp_indexes[i + 1]
            particles_valid = tmp_particles[i] == "E" and tmp_particles[i + 1] == "O"
            particles_valid = particles_valid or (tmp_particles[i] == "O" and tmp_particles[i + 1] == "E")

            # Check if oxygen diffusion can happen.
            if particles_valid and indexes_valid:
                tmp_state_0[i], tmp_state_0[i + 1] = tmp_state_0[i + 1], tmp_state_0[i]
                append_to_final_states(tuple(tmp_state_0))

        return final_states

    # --------------------------------------------------------------------------
    # Process Reaction Methods.
    # --------------------------------------------------------------------------

    def process_reaction_er_reaction(self, initial_state: tuple):
        """ Given a state, returns the list of states that are generated by an
            Elay-Rideal type of reaction. This reaction requires a single site
            with oxygen.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of an Elay-Rideal reaction.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def append_to_final_states(state0: tuple):
            """ Appends the given state to the final state list, with the
                correct format.

                :param state0: The state to be appended.
            """

            final_state = tuple((cp.deepcopy(entry), cp.deepcopy(tmp_indexes[j])) for j, entry in enumerate(state0))
            final_states.append(final_state)

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        final_states = []
        if len(initial_state) < self.o_coo_er:
            return final_states

        # Get the list of states and the numbering.
        tmp_particles, tmp_indexes = self.get_state_elements(initial_state)
        for i, state in enumerate(tmp_particles):
            tmp_state_0 = cp.deepcopy(list(tmp_particles))
            if state == 'O':
                tmp_state_0[i] = "E"
                append_to_final_states(tuple(tmp_state_0))

        return final_states

    def process_reaction_lh_reaction(self, initial_state: tuple):
        """ Given a state, returns the list of states that are generated by a
            Langmuir-Hinshelwood type of reaction. This reaction requires two
            adjacent sites, one with carbon monoxide and the other one with
            oxygen.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of a Langmuir-Hinshelwood reaction.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def append_to_final_states(state0: tuple):
            """ Appends the given state to the final state list, with the
                correct format.

                :param state0: The state to be appended.
            """

            final_state = tuple((cp.deepcopy(entry), cp.deepcopy(tmp_indexes[j])) for j, entry in enumerate(state0))
            final_states.append(final_state)

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        final_states = []
        if len(initial_state) < self.o_coo_lh:
            return final_states

        # Get the list of states and the numbering.
        tmp_particles, tmp_indexes = self.get_state_elements(initial_state)

        # Get the resulting states.
        for i, _ in enumerate(tmp_particles):
            if i == len(tmp_particles) - 1:
                break

            # Always get a copy of the initial state first.
            tmp_state_0 = cp.deepcopy(list(tmp_particles))
            indexes_valid = tmp_indexes[i] + 1 == tmp_indexes[i + 1]
            particles_valid = tmp_particles[i] == "CO" and tmp_particles[i + 1] == "O"
            particles_valid = particles_valid or (tmp_particles[i] == "O" and tmp_particles[i + 1] == "CO")

            # Check if the Langmuir-Hinshelwood reaction can happen.
            if particles_valid and indexes_valid:
                tmp_state_0[i], tmp_state_0[i + 1] = "E", "E"
                append_to_final_states(tuple(tmp_state_0))

        return final_states

    # --------------------------------------------------------------------------
    # Reduce Methods.
    # --------------------------------------------------------------------------

    def reduce_to_unique_states(self, state_dictionary, target_state):
        """ Given a dictionary of states, it attempts to contract all the keys.

            :param state_dictionary: The dictionary of states to be reduced.

            :param target_state: The state that is being targeted to appear in
            the reduced list.

            :return: A list of the reduced states in the format
                    (state, multiplicity).
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def add_remove_states(combinations0):
            """ From the list of UNIQUE states to be examined for contraction,
                returns the index to be added and the indexes to be contracted.
                Contractions are preferred occur in the last index.

                :param: combinations0: The combinations of the states in lists of
                length len(self.states).

                :return: A 2-tuple where the first entry is the state that
                resulted from the contraction and the other entry is the list of
                states that generated the contracted state.
            """

            # Try to contract on the right end.
            contracted_states0 = list(map(lambda x: self.get_contracted_state(x, -1), combinations0))
            contracted_states0 = [contracted for contracted in contracted_states0 if len(contracted[0]) > 0]

            for contracted_state0 in contracted_states0:
                if subindexes_in_subtates(target_state, contracted_state0[0]):
                    return contracted_state0[0], contracted_state0[1]

            # Try to contract on the left end.
            contracted_states0 = list(map(lambda x: self.get_contracted_state(x, 0), combinations0))
            contracted_states0 = [contracted0 for contracted0 in contracted_states0 if len(contracted0[0]) > 0]
            for contracted_state0 in contracted_states0:
                if subindexes_in_subtates(target_state, contracted_state0[0]):
                    return contracted_state0[0], contracted_state0[1]

            # Contraction not possible at either end.
            return [], []

        def subindexes_in_subtates(state0, state1):
            """ Determines if the indexes of state0 are a subset o the indexes
                in state1.

                :param state0: The state whose indexes are to be verified as
                being a subset of the indexes of state 1.

                :param state1: The state whose indexes are to be verified
                against those of state0.

                :return: True, if the indexes of state0 are a proper subset of
                the indexes of state1. False, otherwise.
            """

            indexes0 = set(index[1] for index in state0)
            indexes1 = set(index[1] for index in state1)

            return indexes0.issubset(indexes1)

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        keys = state_dictionary.keys()
        state_dictionary0 = cp.deepcopy(state_dictionary)

        for key in keys:
            length_after, length_before, enter = 0, 0, True

            while not length_after == length_before or enter:
                enter = False
                unique_elements = list(set(state_dictionary0[key]))
                length_before = len(unique_elements)
                combinations = list(itertools.combinations(unique_elements, len(self.states)))
                if len(combinations) == 0:
                    break

                state_to_add, states_to_remove = add_remove_states(combinations)
                if len(state_to_add) > 0:
                    for state_to_remove in states_to_remove:
                        while state_dictionary0[key].count(state_to_remove) > 0:
                            state_dictionary0[key].remove(state_to_remove)

                # Get the unique elements again.
                unique_elements = list(set(state_dictionary0[key]))
                length_after = len(unique_elements)

        return state_dictionary0

    # --------------------------------------------------------------------------
    # Save Methods.
    # --------------------------------------------------------------------------

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
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def save_to_file(save_string0):
            """ Saves the given string to given path.

                :param save_string0: The string to be saved.
            """

            # Auxiliary variables.
            i0 = 0

            # Format the path properly.
            save_path0 = save_path
            save_path0 += "" if save_path[-1] == os.sep else os.sep

            # Get the full save path.
            path0 = save_path0 + file_name + ".txt"

            # Check that the file does not overwite the new file.
            while os.path.isfile(path0):
                # Add a numbered path.
                path0 = save_path0 + file_name + str(i0) + ".txt"

                # Update the counter.
                i0 += 1

            # Open the file to write.
            with open(path0, "w") as fl:
                # Write the string.
                fl.write(save_string0)

        def validate_parameters():
            """ Validates that the parameters are consistent.

                :return save_path0: The valid an correct save path for the file.
                In case there is no save path parameter, it defaults to the
                current folder.
            """

            # Validate the order.
            if not isinstance(order, (int,)):
                raise TypeError("The order parameter must be an integer.")

            # Get a proper file path.
            cwd = os.path.dirname(__file__)
            cwd += "" if cwd[-1] == os.sep else os.sep
            save_path0 = cwd if save_path is None else save_path.strip()

            # Check that it is, indeed, a directory.
            if not os.path.isdir(save_path0):
                raise ValueError(f"The entered path is not valid, i.e., NOT a directory: {save_path0}")

            return save_path0

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        # Validate the parameters before proceeding.
        save_path = validate_parameters()

        # Check if there are equations.
        if len(self.equations) == 0:
            # Give the user a message.
            print("There are currently no equations to show.")
            return

        # Auxiliary variables.
        equation_strings = []
        constraint_strings = []
        initial_conditions_strings = []
        rates_value_strings = []
        raw_state_strings = []
        keys = list(key for key in self.equations[0][1].keys())

        # Get the requested formatter.
        format0 = format_type.strip().lower()
        formatter0 = FormatterManager.get_formatter(format0)

        # For every equation.
        for i, equation in enumerate(self.equations):
            # Get the particular equation.
            equation_strings.append(formatter0.format_equation(equation, order))

            # Get the initial conditions
            initial_conditions_strings.append(formatter0.format_initial_condition(equation[0]))

            # Get the raw state strings.
            raw_state_strings.append(formatter0.format_state(equation[0], raw=True))

        # For every rate.
        for key in keys:
            # Get the rate value.
            rates_value_strings.append(formatter0.format_rate(key, 0.0))

        # For every constraint.
        for i, constraint in enumerate(self.constraints):
            # Get the constraint equations.
            constraint_strings.append(formatter0.format_constraint(constraint))

        # Quantities that will be formatted.
        format_quantities = {
            "constraints": constraint_strings,
            "equations": equation_strings,
            "initial conditions": initial_conditions_strings,
            "rate values": rates_value_strings,
            "raw states": raw_state_strings
        }

        # Get the string to save.
        save_string = formatter0.format_final(format_quantities)

        # Generate the file.
        save_to_file(save_string)

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Constructors and Dunder Methods.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Constructor.
    # --------------------------------------------------------------------------

    def __init__(self, sites: int = 1):
        """ Builds a class that writes the equations for the carbon monoxide -
            oxygen associative reaction on ruthenium (111);
            J. Chem. Phys. 143, 204702 (2015).

            Initializes the class with the standard parameters.

            :param sites: The number of sites that the system has.
        """

        # Initialize the super class.
        super(EquationGenerator, self).__init__(sites, ('CO', 'O', 'E'))

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

        # Strings that represent the rate constants of oxygen related processes.
        self.k_o_ads = "k.O.ads"
        self.k_o_des = "k.O.des",
        self.k_o_dif = "k.O.Dif"

        # Strings that represent the rate constants of carbon monoxide
        # processes.
        self.k_co_ads = "k.CO.ads"
        self.k_co_des = "k.CO.des"
        self.k_co_dif = "k.CO.dif"

        # Strings that represent the rate constants of carbon monoxide - oxygen
        # reaction related processes.
        self.k_coo_lh = "k.COO.lh"
        self.k_coo_er = "k.COO.er"
