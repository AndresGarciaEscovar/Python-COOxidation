""" Contains the tests for the COOxidationKMC class."""

# ------------------------------------------------------------------------------
# Imports.
# ------------------------------------------------------------------------------

# Imports: The unittest module.
import copy
import random
import unittest

from typing import Union
from itertools import product

# Imports: User-defined.
from coOxidation.Program.KMC.COOxidation_KMC import COOxidationKMC
from coOxidation.Program.KMC.COOxidation_parameters import COOxidationKMCParameters


# ------------------------------------------------------------------------------
# Classes.
# ------------------------------------------------------------------------------


class TestCOOxidationKMC(unittest.TestCase):
    """ Tests that the different functions of the COOxidationKMC class are
        working properly.
    """

    # --------------------------------------------------------------------------
    # Auxiliary Methods.
    # --------------------------------------------------------------------------

    @staticmethod
    def get_simulation(length: int = 3, maximum_counter: Union[int, float] = 100) -> tuple:
        """ Returns a fully instantiated and initialized simulation.

            :param length: The length of the lattice.

            :param maximum_counter: The maximum counter of the simulation.

            :return: A fully instantiated and initialized COOxidationKMC
             simulation object.
        """

        parameters = COOxidationKMCParameters()
        parameters.length = length
        parameters.maximum_counter = maximum_counter

        simulation = COOxidationKMC(parameters)

        return simulation, parameters

    def validate_states(self, lattice: list, final: tuple) -> None:
        """ Validates that the final states are consistent with the process.

            :param lattice: The lattice of particles to be validated.

            :param final: The tuple with the final states after undergoing a
             process.
        """

        for i, component in enumerate(lattice):
            self.assertEqual(component, final[i])

    # --------------------------------------------------------------------------
    # Test Methods.
    # --------------------------------------------------------------------------

    def test_choose_process(self):
        """ Tests that the processes are consistently chosen."""

        simulation, parameters = TestCOOxidationKMC.get_simulation()

        for i in range(10_000):
            self.assertTrue(0 <= simulation.choose_process() < len(simulation.rates))

    def test_constructor(self):
        """  Tests that the constructor initially has the required values."""

        simulation, parameters = TestCOOxidationKMC.get_simulation()  # type: COOxidationKMC, COOxidationKMCParameters

        self.assertEqual(len(simulation), parameters.length)
        self.assertEqual(len(simulation.lattice), len(simulation.statistics))

        self.assertEqual(len(simulation.lattice), parameters.length)
        self.assertTrue(all(state == 'E' for state in simulation.lattice))

        self.assertEqual(len(simulation.statistics), parameters.length)
        self.assertTrue(all({'E', 'O', 'CO'} == set(statistics.keys()) for statistics in simulation.statistics))

        for statistics in simulation.statistics:
            for key in statistics.keys():
                self.assertEqual(statistics[key], 0)

        self.assertEqual(simulation.counter_steps, 0)
        self.assertEqual(simulation.counter_time, 0.0)
        self.assertIsInstance(simulation.counter_maximum, (int,))

        self.assertEqual(simulation.seed, parameters.seed)

    def test_get_floats(self):
        """  Tests that the site validation method is working."""

        length = 10
        simulation, parameters = self.get_simulation(length=length)

        for i in range(100_000):
            base = random.choice([123.123612387362, -12.123971237])
            number = simulation.get_float(base)
            self.assertTrue(0 < abs(number) < abs(base))

    def test_get_rates(self):
        """  Tests that the site validation method is working."""

        length = 10
        simulation, parameters = self.get_simulation(length=length)

        for i, rate in enumerate(simulation.rates):
            if i == 0:
                self.assertGreaterEqual(rate, 0.0)
                continue
            self.assertLessEqual(simulation.rates[i - 1], rate)

    def test_get_sites(self):
        """ Tests that the get sites method is working properly."""

        simulation, parameters = self.get_simulation()

        for i in range(10_000):
            sites = simulation.get_sites()
            self.assertEqual(len(sites), 1)
            self.assertTrue(0 <= sites[0] < len(simulation))

        for i in range(10_000):
            sites = simulation.get_sites(two=True)
            self.assertEqual(len(sites), 2)
            self.assertTrue(0 <= sites[0] < len(simulation))
            self.assertTrue(sites[0] == sites[1] + 1 or sites[0] == sites[1] - 1)

    def test_process_adsorb(self):
        """ Tests that the adsorb processes are consistently carried out."""

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def empty_lattice(lattice0: list, state0: tuple = tuple()) -> None:
            """ Sets the lattice to the given state. If the state is empty, the
                lattice will be emptied.

                :param lattice0: The lattice whose state is to be set.

                :param state0: The state to which the respective lattice
                 components must be set. Must match the dimensions.
            """

            if state0 == tuple():
                for i in range(len(lattice0)):
                    lattice0[i] = 'E'
                return

            for i, state0_ in enumerate(state0):
                lattice0[i] = state0_

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        simulation, parameters = TestCOOxidationKMC.get_simulation()
        lattice = simulation.lattice
        states0 = ['CO', 'E', 'O']
        states1 = ['CO', 'O']
        states2 = ['O', 'CO']

        # Desorption of Oxygen.
        sites_list = [[0, 1], [1, 0]]
        for sites in sites_list:
            for particle in states0:
                empty_lattice(lattice, ('E', 'E', particle))
                simulation.process_adsorb(sites, 0)
                self.validate_states(lattice, ('O', 'O', particle))

            for particle in states0:
                for particle1 in states1:
                    empty_lattice(lattice, ('E', particle1, particle))
                    simulation.process_adsorb(sites, 0)
                    self.validate_states(lattice, ('E', particle1, particle))

                for particle1 in states1:
                    empty_lattice(lattice, (particle1, 'E', particle))
                    simulation.process_adsorb(sites, 0)
                    self.validate_states(lattice, (particle1, 'E', particle))

                for particle1, particle2 in product(states1, states1):
                    empty_lattice(lattice, (particle1, particle2, particle))
                    simulation.process_adsorb(sites, 0)
                    self.validate_states(lattice, (particle1, particle2, particle))

        sites_list = [[1, 2], [2, 1]]
        for sites in sites_list:
            for particle in states0:
                empty_lattice(lattice, (particle, 'E', 'E'))
                simulation.process_adsorb(sites, 0)
                self.validate_states(lattice, (particle, 'O', 'O'))

            for particle in states0:
                for particle1 in states1:
                    empty_lattice(lattice, (particle, 'E', particle1))
                    simulation.process_adsorb(sites, 0)
                    self.validate_states(lattice, (particle, 'E', particle1))

                for particle1 in states1:
                    empty_lattice(lattice, (particle, particle1, 'E'))
                    simulation.process_adsorb(sites, 0)
                    self.validate_states(lattice, (particle, particle1, 'E'))

                for particle1, particle2 in product(states1, states1):
                    empty_lattice(lattice, (particle, particle1, particle2))
                    simulation.process_adsorb(sites, 0)
                    self.validate_states(lattice, (particle, particle1, particle2))

        # Desorption of CO.
        sites = [0]
        for particle1, particle2 in product(states0, states0):
            empty_lattice(lattice, ('E', particle1, particle2))
            simulation.process_adsorb(sites, 3)
            self.validate_states(lattice, ('CO', particle1, particle2))
            for particle3 in states2:
                empty_lattice(lattice, (particle3, particle1, particle2))
                simulation.process_adsorb(sites, 3)
                self.validate_states(lattice, (particle3, particle1, particle2))

        sites = [1]
        for particle1, particle2 in product(states0, states0):
            empty_lattice(lattice, (particle1, 'E', particle2))
            simulation.process_adsorb(sites, 3)
            self.validate_states(lattice, (particle1, 'CO', particle2))
            for particle3 in states2:
                empty_lattice(lattice, (particle1, particle3, particle2))
                simulation.process_adsorb(sites, 3)
                self.validate_states(lattice, (particle1, particle3, particle2))

        sites = [2]
        for particle1, particle2 in product(states0, states0):
            empty_lattice(lattice, (particle1, particle2, 'E'))
            simulation.process_adsorb(sites, 3)
            self.validate_states(lattice, (particle1, particle2, 'CO'))
            for particle3 in states2:
                empty_lattice(lattice, (particle1, particle2, particle3))
                simulation.process_adsorb(sites, 3)
                self.validate_states(lattice, (particle1, particle2, particle3))

    def test_process_desorb(self):
        """ Tests that the adsorb processes are consistently carried out."""

        # //////////////////////////////////////////////////////////////////////
        # Auxiliary functions.
        # //////////////////////////////////////////////////////////////////////

        def empty_lattice(lattice0: list, state0: tuple = tuple()) -> None:
            """ Sets the lattice to the given state. If the state is empty, the
                lattice will be emptied.

                :param lattice0: The lattice whose state is to be set.

                :param state0: The state to which the respective lattice
                 components must be set. Must match the dimensions.
            """

            if state0 == tuple():
                for i in range(len(lattice0)):
                    lattice0[i] = 'E'
                return

            for i, state0_ in enumerate(state0):
                lattice0[i] = state0_

        # //////////////////////////////////////////////////////////////////////
        # Implementation.
        # //////////////////////////////////////////////////////////////////////

        simulation, parameters = TestCOOxidationKMC.get_simulation()
        lattice = simulation.lattice
        states0 = ['CO', 'E', 'O']
        states1 = ['CO', 'E']
        states2 = ['O', 'E']

        # Desorption of Oxygen.
        sites_list = [[0, 1], [1, 0]]
        for sites in sites_list:
            for particle in states0:
                empty_lattice(lattice, ('O', 'O', particle))
                simulation.process_desorb(sites, 1)
                self.validate_states(lattice, ('E', 'E', particle))

            for particle in states0:
                for particle1 in states1:
                    empty_lattice(lattice, ('O', particle1, particle))
                    simulation.process_desorb(sites, 1)
                    self.validate_states(lattice, ('O', particle1, particle))

                for particle1 in states1:
                    empty_lattice(lattice, (particle1, 'O', particle))
                    simulation.process_desorb(sites, 1)
                    self.validate_states(lattice, (particle1, 'O', particle))

                for particle1, particle2 in product(states1, states1):
                    empty_lattice(lattice, (particle1, particle2, particle))
                    simulation.process_desorb(sites, 1)
                    self.validate_states(lattice, (particle1, particle2, particle))

        sites_list = [[1, 2], [2, 1]]
        for sites in sites_list:
            for particle in states0:
                empty_lattice(lattice, (particle, 'O', 'O'))
                simulation.process_desorb(sites, 1)
                self.validate_states(lattice, (particle, 'E', 'E'))

            for particle in states0:
                for particle1 in states1:
                    empty_lattice(lattice, (particle, 'O', particle1))
                    simulation.process_desorb(sites, 1)
                    self.validate_states(lattice, (particle, 'O', particle1))

                for particle1 in states1:
                    empty_lattice(lattice, (particle, particle1, 'O'))
                    simulation.process_desorb(sites, 1)
                    self.validate_states(lattice, (particle, particle1, 'O'))

                for particle1, particle2 in product(states1, states1):
                    empty_lattice(lattice, (particle, particle1, particle2))
                    simulation.process_desorb(sites, 1)
                    self.validate_states(lattice, (particle, particle1, particle2))

        # Desorption of CO.
        sites = [0]
        for particle1, particle2 in product(states0, states0):
            empty_lattice(lattice, ('CO', particle1, particle2))
            simulation.process_desorb(sites, 4)
            self.validate_states(lattice, ('E', particle1, particle2))
            for particle3 in states2:
                empty_lattice(lattice, (particle3, particle1, particle2))
                simulation.process_desorb(sites, 4)
                self.validate_states(lattice, (particle3, particle1, particle2))

        sites = [1]
        for particle1, particle2 in product(states0, states0):
            empty_lattice(lattice, (particle1, 'CO', particle2))
            simulation.process_desorb(sites, 4)
            self.validate_states(lattice, (particle1, 'E', particle2))
            for particle3 in states2:
                empty_lattice(lattice, (particle1, particle3, particle2))
                simulation.process_desorb(sites, 4)
                self.validate_states(lattice, (particle1, particle3, particle2))

        sites = [2]
        for particle1, particle2 in product(states0, states0):
            empty_lattice(lattice, (particle1, particle2, 'CO'))
            simulation.process_desorb(sites, 4)
            self.validate_states(lattice, (particle1, particle2, 'E'))
            for particle3 in states2:
                empty_lattice(lattice, (particle1, particle2, particle3))
                simulation.process_desorb(sites, 4)
                self.validate_states(lattice, (particle1, particle2, particle3))

    def test_update_counters(self):
        """  Tests that the site validation method is working."""

        length = 10
        simulation, parameters = self.get_simulation(length=length)

        counter_steps = -1
        for i in range(30):
            counter_steps = copy.deepcopy(simulation.counter_steps)
            counter_time = copy.deepcopy(simulation.counter_time)
            simulation.update_counters()

            self.assertLess(counter_steps, simulation.counter_steps)
            self.assertLess(counter_time, simulation.counter_time)

        self.assertTrue(counter_steps, 29)

    def test_validate_in_lattice(self):
        """  Tests that the site validation method is working."""

        length = 10
        simulation, parameters = self.get_simulation(length=length)

        valid = [i for i in range(len(simulation))]
        for i in valid:
            self.assertTrue(simulation.validate_in_lattice(i))

        not_valid = [i for i in range(len(simulation), 2 * len(simulation))]
        not_valid.extend([i for i in range(-len(simulation), 0)])
        for i in not_valid:
            self.assertFalse(simulation.validate_in_lattice(i))


if __name__ == '__main__':
    # Run the tests.
    unittest.main()
