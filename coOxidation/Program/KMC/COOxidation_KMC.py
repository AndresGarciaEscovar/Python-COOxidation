""" Contains the COOxidationKMC class."""

# ------------------------------------------------------------------------------
# Imports.
# ------------------------------------------------------------------------------

# Imports: General.
import numpy

# Imports: User-defined.
from coOxidation.Program.KMC.COOxidation_parameters import COOxidationKMCParameters

# ------------------------------------------------------------------------------
# Classes.
# ------------------------------------------------------------------------------


class COOxidationKMC:

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
        random, id_ = self.get_float(self.rates[-1]), numpy.inf
        for i, rate in enumerate(self.rates):
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

    @staticmethod
    def get_rates() -> tuple:
        """ Gets the rates of the system. More rates can be added if needed.

            :return: The cumulative rates of the system.
        """

        rates = (
            2 * 1.0,  # Oxygen adsorption rate.
            2 * 1.0,  # Oxygen desorption rate.
            2 * 1.0,  # Oxygen diffusion rate.
            1 * 1.0,  # Carbon monoxide adsorption rate.
            1 * 1.0,  # Carbon monoxide desorption rate.
            2 * 1.0,  # Carbon monoxide diffusion rate.
            2 * 1.0,  # Langmuir-Hinshelwood reaction rate.
            1 * 1.0,  # Elay-Rideal reaction rate.
        )
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

        if not self.lattice[sites[0]] == 'E':
            return

        if process_id == 0:
            if not (self.validate_in_lattice(sites[1]) and self.lattice[sites[1]] == 'E'):
                return

            self.lattice[sites[0]], self.lattice[sites[1]] = 'O', 'O'
            return

        self.lattice[sites[0]] = 'CO'

    def process_desorb(self, sites: list, process_id: int) -> None:
        """ Tries to desorb a carbon monoxide molecule or an oxygen pair.

            :param sites: The list of sites that are involved in the adsorption
             process. The zeroth site will always be a valid site, if there are
             more sites, these must be calculated before hand.

            :param process_id: The id of the process.
        """

        if self.lattice[sites[0]] == 'E':
            return

        if process_id == 1:
            if not self.validate_in_lattice(sites[1]):
                return

            if not (self.lattice[sites[0]] == 'O' and self.lattice[sites[1]] == 'O'):
                return

            self.lattice[sites[0]], self.lattice[sites[1]] = 'E', 'E'
            return

        self.lattice[sites[0]] = 'E' if self.lattice[sites[0]] == 'CO' else self.lattice[sites[0]]

    def process_diffusion(self, sites: list, process_id: int) -> None:
        """ Tries to diffuse a carbon monoxide molecule or an oxygen atom.

            :param sites: The list of sites that are involved in the adsorption
             process. The zeroth site will always be a valid site, if there are
             more sites, these must be calculated before hand.

            :param process_id: The id of the process.
        """

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

        if self.lattice[sites[0]] == 'E':
            return

        if process_id == 7:
            if self.lattice[sites[0]] == 'O':
                self.lattice[sites[0]] = 'E'
            return

        if not self.validate_in_lattice(sites[1]):
            return

        states = ['CO', 'O']
        states.remove(self.lattice[sites[0]])
        if self.lattice[sites[1]] == states[0]:
            self.lattice[sites[0]],  self.lattice[sites[1]] = 'E',  'E'

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

    def run_simulation(self, repetitions: int = 1) -> None:
        """ Run the simulation.

            :param repetitions: The number of times the simulation must be
             run, for statistics purposes.
        """

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def must_continue() -> bool:
            """ Determines if the simulation must continue, i.e., the guard
                variable has surpassed the maximum counter.

                :return:
            """

            if isinstance(self.counter_maximum, (float,)) and self.counter_time <= self.counter_maximum:
                return True

            if isinstance(self.counter_maximum, (int,)) and self.counter_steps <= self.counter_maximum:
                return True

            return False

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        self.reset_simulation(completely=True)

        for _ in range(repetitions):
            self.reset_simulation(completely=False)

            while True:
                self.update_counters()
                if not must_continue():
                    continue

                self.choose_move()
                self.statistics_record()

        self.statistics_record(repetitions=repetitions, normalize=True)

    # --------------------------------------------------------------------------
    # Statistics Methods.
    # --------------------------------------------------------------------------

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
        self.counter_time -= numpy.log(self.get_float(1.0)) / self.rates[-1]

    # --------------------------------------------------------------------------
    # Validate methods.
    # --------------------------------------------------------------------------

    def validate_in_lattice(self, site: int) -> bool:
        """ Validates if the given site is in the lattice.

            :param site: The site to be validated.

            :return: True, if the site is in the lattice. False otherwise.
        """
        return 0 <= site < len(self)

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
        self.rates = COOxidationKMC.get_rates()

        self.counter_steps = 0
        self.counter_time = 0.0
        self.counter_maximum = parameters.maximum_counter

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
