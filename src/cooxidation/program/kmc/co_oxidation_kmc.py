""" Contains the COOxidationKMC class."""

# ------------------------------------------------------------------------------
# Imports.
# ------------------------------------------------------------------------------

# Imports: General.
import copy
import csv
import datetime
import numpy

from datetime import datetime

# Imports: User-defined.
from co2oxidation.program.kmc.cooxidation_parameters import cooxidationkmcparameters

# ------------------------------------------------------------------------------
# Classes.
# ------------------------------------------------------------------------------


class COOxidationKMC:
    """ Class that carries out the simulation.

        Parameters:

        - self.counter_maximum: The maximum time or steps that must be
          simulated in a single simulation.

        - self.counter_steps: The counter that keeps track of the step progress
          of a single simulation.

        - self.counter_time: The counter that keeps track of the time progress
          of a single simulation.

        - self.generator = The specific random number generator from where to
          get random numbers.

        - self.lattice: A list of strings that denotes the state of the system
          at a given instant.

        - self.rates = A dictionary with the rates of the system.

        - self.rates_cumulative: The list that contains the cumulative rates of
          the system, in the specific order given in the
          COOxidationKMCParameters data class.

        - self.repetitions: The number of simulations over which to average.

        - self.seed = The seed with which to seed the random number generator.

        - self.statistics: A list of dictionaries that contains the statistics
          of the simulation.
    """

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Methods.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Choose Methods.
    # --------------------------------------------------------------------------

    def choose_process(self) -> int:
        """ Randomly chooses a process from the available processes.

            :return: The id of the process to be executed.
        """
        random, id_ = self.get_float(self.rates_cumulative[-1]), numpy.inf
        for i, rate in enumerate(self.rates_cumulative):
            if random < rate:
                return i

        raise ValueError("A number within the processes must be returned")

    def choose_move(self) -> None:
        """ Chooses a random move to be performed.
        """

        # Get the site and process
        process_id = self.choose_process()

        # Perform an adsorption/desorption/diffusion/reaction move.
        if process_id in [0, 3]:  # Adsorption.
            two = process_id == 0
            self.process_adsorb(self.get_sites(two=two), process_id)
            return

        if process_id in [1, 4]:  # Desorption.
            two = process_id == 1
            self.process_desorb(self.get_sites(two=two), process_id)
            return

        if process_id in [2, 5]:  # Diffusion.
            self.process_diffusion(self.get_sites(two=True), process_id)
            return

        if process_id in [6, 7]:  # Reaction
            two = process_id == 6
            self.process_reaction(self.get_sites(two=two), process_id)
            return

        raise ValueError(f"The selected rate must be between 0 and 7. Selected process: {process_id}")

    # --------------------------------------------------------------------------
    # Format Methods.
    # --------------------------------------------------------------------------

    def format_columns(self) -> list:
        """ Gets the columns information.

            :return: A list of lists that contains the columns information.
        """

        str0 = [["Site\\Probability", "E", "O", "CO"]]

        statiscs = copy.deepcopy(self.statistics)
        for i in range(len(self)):
            line = [
                f"{i + 1}",
                *list(f"{statiscs[i][key]:.7f}" for key in ['E', 'O', 'CO'])
            ]
            str0.append(line)

        return str0

    def format_header(self) -> list:
        """ Gets the main information of the simulation.

            :return: A list with the information of the simulation.
        """

        time_stamp = datetime.now().strftime("20%y-%m-%d--%H-%M-%S")
        counter_type = "Steps" if isinstance(self.counter_maximum, (int,)) else "Time"
        rates_strings = [f"{key} = {self.rates[key]:f}" for key in self.rates.keys()]

        str0 = [
            "CO Oxidation on Ru(111)",
            f"date--time = {time_stamp}",
            f"Repetitions = {self.repetitions}",
            f"Simulation Steps = {self.counter_steps}",
            f"Simulation Time = {self.counter_time}",
            f"Maximum {counter_type} Counter = {self.counter_maximum}",
            f"Generator Seed = {self.seed}",
            *rates_strings
        ]

        return str0

    # --------------------------------------------------------------------------
    # Get Methods.
    # --------------------------------------------------------------------------

    def get_float(self, base: float) -> float:
        """ Gets a random number in the range (lower_0, upper_0)

            :param  base: The base number to return a random floating point
             number x in the range 0 < x < base. If the base is zero, it returns
             zero.
        """

        if base == 0.0:
            return 0.0

        number = self.generator.random()
        while not (0.0 < number < 1.0):
            number = self.generator.random()

        return number * base

    def get_rates(self) -> tuple:
        """ Gets the rates of the system. More rates can be added if needed.

            :return: The cumulative rates of the system.
        """

        rates = [self.rates[key] for key in self.rates.keys()]
        rates = tuple(sum(rates[0:i + 1]) for i, _ in enumerate(rates))
        return rates

    def get_sites(self, two=False) -> list:
        """ Gets a list with a single site. If a second one is requested, it
            returns the neighbor site towards the left or right, with equal
            probability.

            :param two: If two sites are requested.

            :return: A list of sites.
        """

        sites = [self.generator.integers(0, len(self))]  # type: list
        if two:
            site = sites[0] + self.generator.choice([-1, 1])
            sites.append(site)

        return sites

    # --------------------------------------------------------------------------
    # Process Methods.
    # --------------------------------------------------------------------------

    def process_adsorb(self, sites: list, process_id: int) -> None:
        """ Tries to adsorb a carbon monoxide molecule or an oxygen pair.

            :param sites: The list of sites that are involved in the adsorption
             process. The zeroth site will always be a valid site, if there are
             more sites, these must be calculated before hand.

            :param process_id: The id of the process.
        """

        COOxidationKMC.validate_process(process_id, (0, 3))

        if not self.lattice[sites[0]] == 'E':
            return

        if process_id == 0:
            if not (self.validate_in_lattice(sites[1]) and self.lattice[sites[1]] == 'E'):
                return

            self.lattice[sites[0]], self.lattice[sites[1]] = 'O', 'O'
            return

        else:
            self.lattice[sites[0]] = 'CO'
            return

    def process_desorb(self, sites: list, process_id: int) -> None:
        """ Tries to desorb a carbon monoxide molecule or an oxygen pair.

            :param sites: The list of sites that are involved in the adsorption
             process. The zeroth site will always be a valid site, if there are
             more sites, these must be calculated before hand.

            :param process_id: The id of the process.
        """

        COOxidationKMC.validate_process(process_id, (1, 4))

        if self.lattice[sites[0]] == 'E':
            return

        if process_id == 1:
            if not self.validate_in_lattice(sites[1]):
                return

            if not (self.lattice[sites[0]] == 'O' and self.lattice[sites[1]] == 'O'):
                return

            self.lattice[sites[0]], self.lattice[sites[1]] = 'E', 'E'
            return

        else:
            self.lattice[sites[0]] = 'E' if self.lattice[sites[0]] == 'CO' else self.lattice[sites[0]]
            return

    def process_diffusion(self, sites: list, process_id: int) -> None:
        """ Tries to diffuse a carbon monoxide molecule or an oxygen atom.

            :param sites: The list of sites that are involved in the adsorption
             process. The zeroth site will always be a valid site, if there are
             more sites, these must be calculated before hand.

            :param process_id: The id of the process.
        """

        COOxidationKMC.validate_process(process_id, (2, 5))

        if self.lattice[sites[0]] == 'E':
            return

        if not (self.validate_in_lattice(sites[1]) and self.lattice[sites[1]] == 'E'):
            return

        if (process_id == 2 and self.lattice[sites[0]] == 'O') or (process_id == 5 and self.lattice[sites[0]] == 'CO'):
            self.lattice[sites[0]],  self.lattice[sites[1]] = self.lattice[sites[1]],  self.lattice[sites[0]]

    def process_reaction(self, sites: list, process_id: int) -> None:
        """ Tries to perform a dissociative reaction between carbon monoxide
            and oxygen.

            :param sites: The list of sites that are involved in the adsorption
             process. The zeroth site will always be a valid site, if there are
             more sites, these must be calculated before hand.

            :param process_id: The id of the process.
        """

        COOxidationKMC.validate_process(process_id, (6, 7))

        if self.lattice[sites[0]] == 'E':
            return

        if process_id == 7:
            if self.lattice[sites[0]] == 'O':
                self.lattice[sites[0]] = 'E'
            return

        else:
            if not self.validate_in_lattice(sites[1]):
                return

            states = ['CO', 'O']
            states.remove(self.lattice[sites[0]])
            if self.lattice[sites[1]] == states[0]:
                self.lattice[sites[0]],  self.lattice[sites[1]] = 'E',  'E'

            return

    # --------------------------------------------------------------------------
    # Reset Methods.
    # --------------------------------------------------------------------------

    def reset_simulation(self, completely=False) -> None:
        """ Resets the simulation parameters to their initial value.

            :param completely: If the simulation must be reset completely.
        """

        self.counter_steps = 0
        self.counter_time = 0.0
        self.lattice = ["E" for _ in range(len(self))]

        if completely:
            self.generator = numpy.random.default_rng(self.seed)
            self.statistics = [{"E": 0, "O": 0, "CO": 0} for _ in range(len(self))]

    # --------------------------------------------------------------------------
    # Run Methods.
    # --------------------------------------------------------------------------

    def run_simulation(self) -> None:
        """ Run the simulation."""

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def must_continue() -> bool:
            """ Determines if the simulation must continue, i.e., the guard
                variable has surpassed the maximum counter.

                :return:
            """

            if isinstance(self.counter_maximum, (float,)) and self.counter_time < self.counter_maximum:
                return True

            if isinstance(self.counter_maximum, (int,)) and self.counter_steps <= self.counter_maximum:
                return True

            return False

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        self.reset_simulation(completely=True)

        for _ in range(self.repetitions):
            self.reset_simulation(completely=False)

            while True:
                self.update_counters()
                if not must_continue():
                    break

                self.choose_move()
            self.statistics_record()

        self.statistics_record(repetitions=self.repetitions, normalize=True)

    # --------------------------------------------------------------------------
    # Statistics Methods.
    # --------------------------------------------------------------------------

    def statistics_save(self, file_name: str, mode: str = 'w') -> None:
        """ Saves the final statistics to the given file.

            :param file_name: The name of the file where the results must be
             saved.

            :param mode: The saving mode; i.e., write, append, etc.
        """

        str_ = [self.format_header(), *self.format_columns(), ["-" * 80]]

        with open(file_name, mode, newline="\n") as fl:
            writer = csv.writer(fl, delimiter=",")
            for line in str_:
                writer.writerow(line)

    def statistics_record(self, repetitions: int = 1, normalize: bool = False) -> None:
        """ Records the statistics for the process.

            :param repetitions: The number of times the simulation has been run.

            :param normalize: True, if the statistics must be normalized. False,
             otherwise.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def normalize_statistics(repetitions0: int) -> None:
            """ Normalizes the quantities for statistical purposes.

                :param repetitions0: The number of times the simulation has been
                 run.
            """
            for i0 in range(len(self)):
                for key0 in ['CO', 'O', 'E']:
                    self.statistics[i0][key0] /= repetitions0

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        if normalize:
            normalize_statistics(repetitions)
            return

        for i in range(len(self)):
            particle = self.lattice[i]
            self.statistics[i][particle] += 1

    # --------------------------------------------------------------------------
    # Update Methods.
    # --------------------------------------------------------------------------

    def update_counters(self) -> None:
        """ Updates the steps and time counters.
        """""
        self.counter_steps += 1
        self.counter_time -= numpy.log(self.get_float(1.0)) / self.rates_cumulative[-1]

    # --------------------------------------------------------------------------
    # Validate methods.
    # --------------------------------------------------------------------------

    def validate_in_lattice(self, site: int) -> bool:
        """ Validates if the given site is in the lattice.

            :param site: The site to be validated.

            :return: True, if the site is in the lattice. False otherwise.
        """
        return 0 <= site < len(self)

    @staticmethod
    def validate_process(process_id: int, processes: tuple) -> None:
        """ Validates that the given process is in the list of available
            processes.

            :param process_id: The process identifier to be validated.

            :param processes: The list of valid process identifiers.
        """

        if process_id not in processes:
            raise ValueError(
                f"The desired process is not valid. Process: {process_id}, "
                f"available processes: {processes}."
            )

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Constructor and Dunder Methods.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Constructor.
    # --------------------------------------------------------------------------

    def __init__(self, parameters: COOxidationKMCParameters):
        """ Constructs a kinetic Monte Carlo simulation for CO oxidation. All
            variables are initially set to a default value.

            :param parameters: The parameters of the simulation.
        """

        self.lattice = ["E" for _ in range(parameters.length)]
        self.statistics = [{"E": 0, "O": 0, "CO": 0} for _ in range(len(self))]
        self.rates = parameters.rates
        self.rates_cumulative = self.get_rates()

        self.counter_steps = 0
        self.counter_time = 0.0
        self.counter_maximum = parameters.maximum_counter
        self.repetitions = parameters.repetitions

        self.seed = parameters.seed
        self.generator = numpy.random.default_rng(self.seed)

    # --------------------------------------------------------------------------
    # Dunder Methods.
    # --------------------------------------------------------------------------

    def __len__(self):
        """ Returns the length of the lattice.

            :return: The length of the lattice.
        """
        return len(self.lattice)
