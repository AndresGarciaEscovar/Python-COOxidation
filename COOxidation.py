# Imports: General.
import copy as cp
import os

import numpy as np
import random


class COOxidationKMC:

    # The current directory where the file is located.
    CURRENT_DIRECTORY = os.path.dirname(__file__)

    # --------------------------------------------------------------------------
    # Getters and Setters.
    # --------------------------------------------------------------------------

    @property
    def configuration_list(self):
        """ Returns the configuration list.

            :return self.__configuration_list: The state of the cells.
        """
        return self.__configuration_list

    @configuration_list.setter
    def configuration_list(self, configuration_list):
        """ Sets the new configuration list.

            :param configuration_list: The new configuration list.
        """
        self.__configuration_list = configuration_list

    @configuration_list.deleter
    def configuration_list(self):
        """ Cannot delete this variable.
        """
        raise AttributeError("Cannot delete this variable, i.e., configuration_list.")

    # --------------------------------------------------------------------------

    @property
    def elapsed_time(self):
        """ Returns the elapsed time.

            :return self.__elapsed time: The elapsed time of the simulation; in
            simulation units.
        """
        return self.__elapsed_time

    @elapsed_time.setter
    def elapsed_time(self, elapsed_time):
        """ Sets the new elapsed time.
            
            :param elapsed_time: The new elapsed time.
        """
        self.__elapsed_time = elapsed_time

    @elapsed_time.deleter
    def elapsed_time(self):
        """ Cannot delete this variable.
        """
        raise AttributeError("Cannot delete this variable, i.e., elapsed_time.")

    # --------------------------------------------------------------------------

    @property
    def final_state_count(self):
        """ Returns the final state count after a simulation.

            :return self.__final_state_count: The list that contains the
            statistics of the simulation.
        """
        return self.__final_state_count

    @final_state_count.setter
    def final_state_count(self, final_state_count):
        """ Sets the new final_state_count.

            :param final_state_count: The list of dictionaries that contains the
            statistics.
        """

        self.__final_state_count = final_state_count

    @final_state_count.deleter
    def final_state_count(self):
        """ Cannot delete this variable.
        """
        raise AttributeError("Cannot delete this variable, i.e., final_state_count.")

    # --------------------------------------------------------------------------

    @property
    def maximum_time(self):
        """ Returns the maximum time the simulation should run for; in
            simulation time units.

            :return self.__maximum_time: The maximum time the simulation should
            run for; in simulation time units.
        """
        return cp.deepcopy(self.__maximum_time)

    @maximum_time.setter
    def maximum_time(self, maximum_time):
        """ Sets the new maximum time.

            :param maximum_time: The new maximum time the simulation should
            run for; in simulation time units.
        """

        try:
            self.__maximum_time

        except AttributeError:
            self.__maximum_time = np.abs(maximum_time)

    @maximum_time.deleter
    def maximum_time(self):
        """ Cannot delete this variable.
        """
        raise AttributeError("Cannot delete this variable, i.e., maximum_time.")

    # --------------------------------------------------------------------------

    @property
    def steps_number(self):
        """ Returns the steps number of the simulation.

            :return self.__steps_number: The list that contains the
            statistics of the simulation.
        """
        return cp.deepcopy(self.__steps_number)

    @steps_number.setter
    def steps_number(self, steps_number):
        """ Sets the new final_state_count.

            :param steps_number: The list of dictionaries that contains the
            statistics.
        """

        try:
            self.__steps_number

        except AttributeError:
            self.__steps_number = np.abs(steps_number)

    @steps_number.deleter
    def steps_number(self):
        """ Cannot delete this variable.
        """
        raise AttributeError("Cannot delete this variable, i.e., steps_number.")

    # ------------------------------------------------------------------------------
    # Get Functions.
    # ------------------------------------------------------------------------------

    def get_rates(self):
        """ Gets the rates of the system. More rates can be added if needed.

            :return rates_0: Returns the cumulative rates of the system.
        """

        # There are 9 rates in total.
        rates_0 = [0 for _ in range(9)]

        # --------------------------------------------------------------------------
        # Oxygen related rates.
        # --------------------------------------------------------------------------

        # Oxygen adsorption, always try to adsorb on neighboring sites (left or
        # right) with a rate of 1.
        rates_0[0] = 2 * 1

        # Oxygen desorption, always try to desorb from neighboring sites (left or
        # right) with a rate of 1.
        rates_0[1] = 2 * 1

        # Oxygen diffusion, two possible moves (left or right) with a rate of 1.
        rates_0[2] = 2 * 1

        # --------------------------------------------------------------------------
        # Carbon monoxide related rates.
        # --------------------------------------------------------------------------

        # Carbon monoxide adsorption.
        rates_0[3] = 1

        # Carbon monoxide desorption.
        rates_0[4] = 1

        # Carbon monoxide diffusion, two possible moves (left or right) with a rate of 1.
        rates_0[5] = 2 * 1

        # --------------------------------------------------------------------------
        # Carbon monoxide - oxygen interaction related rates.
        # --------------------------------------------------------------------------

        # Carbon monoxide - oxygen neighbor reaction on surface, two possible
        # reaction sites (left or right) with a rate of 1.
        rates_0[6] = 2 * 1

        # Carbon monoxide - oxygen reaction in gas (i.e., oxygen on surface binds to
        # CO in gas)
        rates_0[7] = 1

        # Carbon monoxide - oxygen reaction anywhere in lattice, two possible
        # reaction sites (i.e., not the chosen site) with a rate of 1.
        rates_0[8] = 2 * 1

        # --------------------------------------------------------------------------
        # Calculate the cumulative rates.
        # --------------------------------------------------------------------------

        rates_0 = np.array([sum(rates_0[0:i + 1]) for i, _ in enumerate(rates_0)], dtype=np.double)

        return rates_0

    # ------------------------------------------------------------------------------
    # System move functions.
    # ------------------------------------------------------------------------------

    def adsorb_desorb_particle(self, site_0, rate_id_0):
        """ Tries to adsorb a carbon monoxide atom or an oxygen pair, or desorbs a
            an oxygen atom or a carbon monoxide atom.

            :param site_0: The site at which the process will take place.

            :param rate_id_0: The rate that was chosen; determines the process
            adsorption or desorption and the particle type that will undergo the
            process.
        """

        # Check for adsortion/desorption of oxygen.
        if rate_id_0 == 0 or rate_id_0 == 1:

            # Choose a neighbor site.
            site_1 = site_0 + 1 if self.choose_random_number(0, 1) <= 0.5 else site_0 - 1

            # Check that the site is NOT out of bounds.
            if site_1 < 0 or site_1 > 2:
                return

            if rate_id_0 == 0 and self.configuration_list[site_0] == "E" and self.configuration_list[site_1] == "E":
                self.configuration_list[site_0] = "O"
                self.configuration_list[site_1] = "O"
                return

            if rate_id_0 == 1 and self.configuration_list[site_0] == "O" and self.configuration_list[site_1] == "O":
                self.configuration_list[site_0] = "E"
                self.configuration_list[site_1] = "E"
                return

        # Check for adsortion/desorption of carbon monoxide.
        elif rate_id_0 == 3 or rate_id_0 == 4:

            if rate_id_0 == 3 and self.configuration_list[site_0] == "E":
                self.configuration_list[site_0] = "CO"
                return

            if rate_id_0 == 4 and self.configuration_list[site_0] == "CO":
                self.configuration_list[site_0] = "E"
                return

    def move_particle(self, site_0, rate_id_0):
        """ Tries to hop a carbon monoxide atom, or an oxygen atom, to a nearest
            neighbor side with equal probability.

            :param site_0: The site at which the process will take place.

            :param rate_id_0: The rate that was chosen; determines which atom is the
            one that will be hopping.
        """

        # Choose a random direction.
        site_1 = site_0 + 1 if self.choose_random_number(0, 1) <= 0.5 else site_0 - 1

        # Check that the site is NOT out of bounds.
        if site_1 < 0 or site_1 > 2:
            return

        # Conditions under which a swap is possible.
        cond_0 = rate_id_0 == 2 and self.configuration_list[site_0] == "O" and self.configuration_list[site_1] == "E"
        cond_1 = rate_id_0 == 5 and self.configuration_list[site_0] == "CO" and self.configuration_list[site_1] == "E"

        # Make the particle exchange.
        if cond_0 or cond_1:
            tmp_part_0 = cp.deepcopy(self.configuration_list[site_0])
            self.configuration_list[site_0] = self.configuration_list[site_1]
            self.configuration_list[site_1] = tmp_part_0

    def react_particle(self, site_0, rate_id_0):
        """ Reacts a particle or a pair of particles according to the chosen rate,
            i.e., empties the site(s).

            :param site_0: The site at which the process will take place.

            :param rate_id_0: The rate that was chosen; determines which desorption
            process will be attempted.
        """

        # --------------------------------------------------------------------------
        # Auxiliary Functions.
        # --------------------------------------------------------------------------

        def validate_sites(site_1_0, rate_id_1_0):
            """ Validates that the site and reaction index are in the proper range.

                :param site_1_0: The site at which the process will take place.

                :param rate_id_1_0: The rate that was chosen; determines which desorption
                process will be attempted.
            """

            # Validate the new site.
            if site_1_0 not in [0, 1, 2]:
                raise ValueError(f"The requested site for reaction is not valid, it must have a value between 0 and 2."
                                 f"site_0: {site_1_0}.")

            # Validate the reaction index.
            if rate_id_1_0 not in [6, 7, 8]:
                raise ValueError(f"The requested reaction index is not valid, it must have a value between 6 and 8."
                                 f"rate_id_0: {rate_id_1_0}.")

        # --------------------------------------------------------------------------
        # Implementation.
        # --------------------------------------------------------------------------

        # Validate the indexes.
        validate_sites(site_0, rate_id_0)

        # Carbon monoxide - oxygen neighbor reaction on surface, two possible
        # reaction sites (left or right) with a rate of 1.
        if rate_id_0 == 6:
            # Choose a neighboring site.
            site_1 = site_0 + 1 if self.choose_random_number(0, 1) <= 0.5 else site_0 - 1

            # Check that the site is NOT out of bounds.
            if site_1 < 0 or site_1 > 2:
                return

            # Conditions for desorption.
            cond_0 = self.configuration_list[site_0] == "CO" and self.configuration_list[site_1] == "O"
            cond_0 = cond_0 or self.configuration_list[site_0] == "O" and self.configuration_list[site_1] == "CO"

            # Desorb if needed.
            if cond_0:
                self.configuration_list[site_0] = "E"
                self.configuration_list[site_1] = "E"

            return

        # Carbon monoxide - oxygen reaction in gas (i.e., oxygen on surface binds to
        # CO in gas)
        if rate_id_0 == 7:
            # Desorb the oxygen.
            if self.configuration_list[site_0] == "O":
                self.configuration_list[site_0] = "E"

            return

        # Carbon monoxide - oxygen reaction anywhere in lattice, two possible
        # reaction sites (i.e., not the chosen site) with a rate of 1.
        if rate_id_0 == 8:
            # Select the two possible sites from where to choose.
            sites_0 = [i for i in range(3) if not i == site_0]

            # Pick a random site.
            site_1 = random.choice(sites_0)

            # Validate the new site.
            if site_1 < 0 or site_1 > 2 or site_1 == site_0:
                raise ValueError(f"The extra site must be 0, 1 or 2 and cannot be the same as site0."
                                 f"site_0: {site_0}, site_1: {site_1}")

            # Conditions for desorption.
            cond_0 = self.configuration_list[site_0] == "CO" and self.configuration_list[site_1] == "O"
            cond_0 = cond_0 or self.configuration_list[site_0] == "O" and self.configuration_list[site_1] == "CO"

            # Desorb if possible.
            if cond_0:
                self.configuration_list[site_0] = "E"
                self.configuration_list[site_1] = "E"

            return

    # ------------------------------------------------------------------------------
    # Random Functions.
    # ------------------------------------------------------------------------------

    def choose_random_site(self):
        """
            Chooses an integer between 0 and 2.

            :return: An integer number between 0 and 2.
        """
        import random

        return random.choice([int(x) for x in range(3)])

    def choose_random_number(self, lower_0, upper_0):
        """ Gets a random number in the range (lower_0, upper_0)

            :param lower_0: The lower number in the range.

            :param upper_0: The upper number in the range.

            :return random_number: A random number in the range (lower_0, upper_0)
        """
        import random

        # Choose a random number that does NOT include either of the ends.
        range_list = [lower_0, upper_0]
        random_number = random.uniform(min(range_list), max(range_list))
        while random_number == lower_0 or random_number == upper_0:
            random_number = random.uniform(lower_0, upper_0)

        return random_number

    def choose_random_move(self):
        """ Chooses a random move to be performed.
        """

        # Choose a random site.
        site_0 = self.choose_random_site()

        # Choose a random move, i.e., a random rate from the list.
        rates_0 = self.get_rates()
        random_number_0 = self.choose_random_number(0, rates_0[-1])

        # Check that everything is consistent.
        if random_number_0 > rates_0[-1]:
            raise ValueError(f"The number chosen is greater than that of the total rate; this should not happen.\n"
                             f"Maximum rate: {rates_0[-1]}, Chosen rate: {random_number_0}.")

        # Choose the move.
        rate_id_0 = np.inf
        for i_0, rate_0 in enumerate(rates_0):
            if random_number_0 < rate_0:
                rate_id_0 = i_0
                break

        # --------------------------------------------------------------------------
        # Perform an adsorption/desorption move.
        # --------------------------------------------------------------------------

        # Choose the move.
        if rate_id_0 == 0 or rate_id_0 == 3:  # Adsorption.
            self.adsorb_desorb_particle(site_0, rate_id_0)

        elif rate_id_0 == 1 or rate_id_0 == 4:  # Desorption.
            self.adsorb_desorb_particle(site_0, rate_id_0)

        elif rate_id_0 == 2 or rate_id_0 == 5:  # Diffusion.
            self.move_particle(site_0, rate_id_0)

        elif rate_id_0 == 6:  # CO-O reaction on surface.
            self.react_particle(site_0, rate_id_0)

        elif rate_id_0 == 7:  # CO-O reaction in gas.
            self.react_particle(site_0, rate_id_0)

        elif rate_id_0 == 8:  # CO-O reaction anywhere.
            self.react_particle(site_0, rate_id_0)

        else:
            raise ValueError(f"The chosen rate must be between 0 and 8. Current rate_id: {rate_id_0}")

    def choose_time_increment(self):
        """ Increases the elapsed time by a random amount.
        """
        self.elapsed_time = self.elapsed_time - np.log(self.choose_random_number(0, 1)) / (3 * self.get_rates()[-1])

    # --------------------------------------------------------------------------
    # Run a simulation.
    # --------------------------------------------------------------------------

    def run_simulation(self, details=False):
        """ Run the simulation.

            :param details: If details of the simulation must be printed.
        """

        # Do it for the given number of steps.
        for i in range(self.steps_number):

            # Print the details if needed.
            if details:
                print(f"Running simulation {i + 1} of {self.steps_number}")

            # Reset the elapsed time.
            self.elapsed_time = 0.0

            # Reset the state list.
            self.configuration_list = ["E" for _ in range(3)]

            # Time must increment at least once before the simulation begins.
            self.choose_time_increment()

            # While the elapsed time is less than the maximum time.
            while self.elapsed_time < self.maximum_time:
                # Perform a random move.
                self.choose_random_move()

                # Increment the time.
                self.choose_time_increment()

            # Take the statistics.
            for i in range(3):
                self.final_state_count[i][self.configuration_list[i]] += 1

        # Take the statistics.
        for i in range(3):
            for particle in ["E", "O", "CO"]:
                self.final_state_count[i][particle] /= self.steps_number

    # --------------------------------------------------------------------------
    # Save state functions.
    # --------------------------------------------------------------------------

    def save_to_file(self, name="results", append_to_file=True):
        """ Saves the resulst to a text file in the current directory.

            :param name: The name of the file without the extension. That by
            default is results.

            :param append_to_file: If the results should be appended to the
            file instead of written to a different file.
        """

        # Set the base name for the file.
        base_name0 = COOxidationKMC.CURRENT_DIRECTORY + os.sep + name + str(0) + ".txt"

        # Don't overwrite the file.
        i = 0
        while os.path.isfile(base_name0) and not append_to_file:
            base_name0 = COOxidationKMC.CURRENT_DIRECTORY + os.sep + name + str(i) + ".txt"
            i += 1

        # Set the correct mode
        writing_mode = "a" if append_to_file else "w"

        # Save or append to the file.
        with open(base_name0, writing_mode) as fl:

            # Simulation data.
            sim_data = f"Steps: {self.steps_number}, Maximum Time: {self.maximum_time}\n"
            fl.writelines(sim_data)

            # Make a table (to read as a Pandas data frame).
            columns = ["Site", "P(E)", "P(O)", "P(CO)"]
            fl.write("\t".join(columns) + "\n")

            # Gather the results and print them.
            for i in range(3):
                sim_data = [f"{i + 1: 2d}"]
                sim_data += [f"{self.final_state_count[i][j]: .7f}" for j in ["E", "O", "CO"]]
                fl.write("\t".join(sim_data) + "\n")

    # --------------------------------------------------------------------------
    # Constructor.
    # --------------------------------------------------------------------------

    def __init__(self, steps_number, maximum_time):
        """ Constructs a kinetic Monte Carlo simulation for CO oxidation. All
            variables are initially set to a default value.
        """

        # The configuration array will always be set to the standard.
        self.configuration_list = ["E", "E", "E"]

        # Always set the elapsed time to zero.
        self.elapsed_time = np.double(0.0)

        # A dictionary with the particle types and their final count.
        self.final_state_count = [{"E": 0, "O": 0, "CO": 0} for _ in range(3)]

        # Set the number of steps.
        self.steps_number = steps_number

        # Maximum time, i.e., time to run the simulation.
        self.maximum_time = maximum_time


if __name__ == "__main__":
    """ Explanation: To determine the probability of finding a specific particle
        at a site [1 - 3] after a given time t (provided by the user), where t
        is given in units of 1/k; k is the total rate of the system.
        
        To determine the probability, an ensemble average is taken over a given
        number of simulations (provided by the user). Notice that to get good
        statistics we must average over a large number of simulations, with the
        caveat that the more simulations we make to get this average, the more
        time it will take, i.e., we must sacrifice either time or accuracy, or
        find a satisfactory balance.
    """

    """ Instructions:
            * Rates for the different processes must be provided.
                - In the code, look for the get_rates(self) function.
                
                - In the place where the rates are specified, change the
                ones (1) by the corresponding rates. DO NOT touch the
                pre-factors, i.e, the twos (2).
            
            *
    """

    # Create a simulation with parameters:
    #   * n_steps: The number of steps from which to take statistics.
    #
    #   * maximum_simulation_time: The time, in simulation units, for which the
    #   system runs.

    n_steps = 10_000
    maximum_simulation_time = 89

    simulation = COOxidationKMC(10_000, 89)

    # Run the simulation:
    #   * simulation_details: True, if details about the progress of the
    #   simulation should be provided, otherwise False.
    simulation_details = False

    simulation.run_simulation(simulation_details)

    # Saves the data gathered by the simulation:
    #   * append_to_file: True, if the data is to be appended to an existing
    #   file; if the file does not exist, a new one will be created. False, if
    #   a new file should be created for the current run; in the case the file
    #   already exists, it will create a numbered version of the file.
    #
    #   * results_file_name: Name of the file that will store the results. You
    #   can choose not to provide one (delete the parameter from the function)
    #   and it will use the default.
    file_name = "results_file"
    append_to_file_0 = True

    simulation.save_to_file(name=file_name, append_to_file=append_to_file_0)
