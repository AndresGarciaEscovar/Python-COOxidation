"""
    Contains the tests for the COOxidationKMC class.
"""


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Imports.
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# Standard library.
import copy
import random
import unittest

from typing import Union
from itertools import product

# Users.
from coOxidation.Program.KMC.COOxidation_KMC import COOxidationKMC
from coOxidation.Program.KMC.COOxidation_parameters import COOxidationKMCParameters


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Classes.
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


class TestCOOxidationKMC(unittest.TestCase):
    """
        Tests that the different functions of the COOxidationKMC class are
        working properly.
    """
    # /////////////////////////////////////////////////////////////////////////
    # Methods - Auxiliary
    # /////////////////////////////////////////////////////////////////////////

    @staticmethod
    def get_simulation(
        length: int = 3, maximum_counter: Union[int, float] = 100
    ) -> tuple:
        """
            Returns a fully instantiated and initialized simulation.

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
        """
            Validates that the final states are consistent with the process.

            :param lattice: The lattice of particles to be validated.

            :param final: The tuple with the final states after undergoing a
             process.
        """
        for i, component in enumerate(lattice):
            self.assertEqual(component, final[i])

    # /////////////////////////////////////////////////////////////////////////
    # Test Methods.
    # /////////////////////////////////////////////////////////////////////////

    def test_choose_process(self):
        """
            Tests that the processes are consistently chosen.
        """
        simulation, _ = TestCOOxidationKMC.get_simulation()

        for _ in range(10_000):
            process = simulation.choose_process()
            self.assertTrue(0 <= process < len(simulation.rates))

    def test_constructor(self):
        """
            Tests that the constructor initially has the required values.
        """
        simulation, parameters = TestCOOxidationKMC.get_simulation()

        self.assertEqual(len(simulation), parameters.length)
        self.assertEqual(len(simulation.lattice), len(simulation.statistics))

        self.assertEqual(len(simulation.lattice), parameters.length)
        self.assertTrue(all(state == 'E' for state in simulation.lattice))

        self.assertEqual(len(simulation.statistics), parameters.length)
        self.assertTrue(
            all(
                {'E', 'O', 'CO'} == set(statistics.keys())
                for statistics in simulation.statistics
            )
        )

        for statistics in simulation.statistics:
            for key in statistics.keys():
                self.assertEqual(statistics[key], 0)

        self.assertEqual(simulation.counter_steps, 0)
        self.assertEqual(simulation.counter_time, 0.0)
        self.assertIsInstance(simulation.counter_maximum, (int,))

        self.assertEqual(simulation.seed, parameters.seed)

    def test_get_floats(self):
        """
            Tests that the site validation method is working.
        """
        length = 10
        simulation, _ = self.get_simulation(length=length)

        for _ in range(100_000):
            base = random.choice([123.123612387362, -12.123971237])
            number = simulation.get_float(base)
            self.assertTrue(0 < abs(number) < abs(base))

    def test_get_rates(self):
        """
            Tests that the site validation method is working.
        """
        length = 10
        simulation, _ = self.get_simulation(length=length)

        for i, rate in enumerate(simulation.rates_cumulative):
            if i == 0:
                self.assertGreaterEqual(rate, 0.0)
                continue
            self.assertLessEqual(simulation.rates_cumulative[i - 1], rate)

    def test_get_sites(self):
        """
            Tests that the get sites method is working properly.
        """
        simulation, _ = self.get_simulation()

        for _ in range(10_000):
            sites = simulation.get_sites()
            self.assertEqual(len(sites), 1)
            self.assertTrue(0 <= sites[0] < len(simulation))

        for _ in range(10_000):
            sites = simulation.get_sites(two=True)
            self.assertEqual(len(sites), 2)
            self.assertTrue(0 <= sites[0] < len(simulation))
            self.assertTrue(
                sites[0] == sites[1] + 1 or sites[0] == sites[1] - 1
            )

    def test_process_adsorb(self):
        """
            Tests that the adsorb processes are consistently carried out.
        """
        # ---------------------------------------------------------------------
        # Auxiliary functions.
        # ---------------------------------------------------------------------

        def empty_lattice(lattice0: list, state0: tuple = tuple()) -> None:
            """
                Sets the lattice to the given state. If the state is empty, the
                lattice will be emptied.

                :param lattice0: The lattice whose state is to be set.

                :param state0: The state to which the respective lattice
                 components must be set. Must match the dimensions.
            """
            if state0 == tuple():
                for i, _ in enumerate(lattice0):
                    lattice0[i] = 'E'
                return

            for i, state0_ in enumerate(state0):
                lattice0[i] = state0_

        # ---------------------------------------------------------------------
        # Implementation
        # ---------------------------------------------------------------------

        simulation, _ = TestCOOxidationKMC.get_simulation()
        lattice = simulation.lattice
        base_states = ['CO', 'O', 'E']
        states = list(product(*([base_states] * 3)))

        # Oxygen adsorption.
        states_ = [('E', 'E', particle) for particle in base_states]

        for site in [[1, 0], [0, 1]]:
            for state in states:
                empty_lattice(lattice, state)
                simulation.process_adsorb(site, 0)
                if state in states_:
                    self.validate_states(lattice, ('O', 'O', state[2]))
                else:
                    self.validate_states(lattice, state)

        states_ = [(particle, 'E', 'E') for particle in base_states]

        for site in [[1, 2], [2, 1]]:
            for state in states:
                empty_lattice(lattice, state)
                simulation.process_adsorb(site, 0)
                if state in states_:
                    self.validate_states(lattice, (state[0], 'O', 'O'))
                else:
                    self.validate_states(lattice, state)

        # Carbon monoxide adsorption.
        states_ = [
            ('E', *particle)
            for particle in product(*([base_states] * 2))
        ]
        site = [0]
        for state in states:
            empty_lattice(lattice, state)
            simulation.process_adsorb(site, 3)
            if state in states_:
                self.validate_states(lattice, ('CO', *state[1:]))
            else:
                self.validate_states(lattice, state)

        states_ = [
            (particle[0], 'E', particle[1])
            for particle in product(*([base_states] * 2))
        ]
        site = [1]

        for state in states:
            empty_lattice(lattice, state)
            simulation.process_adsorb(site, 3)
            if state in states_:
                self.validate_states(lattice, (state[0], 'CO', state[2]))
            else:
                self.validate_states(lattice, state)

        states_ = [
            (*particle, 'E') for particle in product(*([base_states] * 2))
        ]
        site = [2]

        for state in states:
            empty_lattice(lattice, state)
            simulation.process_adsorb(site, 3)

            if state in states_:
                self.validate_states(lattice, (*state[:-1], 'CO'))
            else:
                self.validate_states(lattice, state)

    def test_process_desorb(self):
        """
            Tests that the desorb processes are consistently carried out.
        """
        # ---------------------------------------------------------------------
        # Auxiliary functions.
        # ---------------------------------------------------------------------

        def empty_lattice(lattice0: list, state0: tuple = tuple()) -> None:
            """ Sets the lattice to the given state. If the state is empty, the
                lattice will be emptied.

                :param lattice0: The lattice whose state is to be set.

                :param state0: The state to which the respective lattice
                 components must be set. Must match the dimensions.
            """
            if state0 == tuple():
                for i, _ in enumerate(lattice0):
                    lattice0[i] = 'E'
                return

            for i, state0_ in enumerate(state0):
                lattice0[i] = state0_

        # ---------------------------------------------------------------------
        # Implementation
        # ---------------------------------------------------------------------

        simulation, _ = TestCOOxidationKMC.get_simulation()
        lattice = simulation.lattice
        base_states = ['CO', 'O', 'E']
        states = list(product(*([base_states] * 3)))

        # Oxygen desorption.
        states_ = [('O', 'O', particle) for particle in base_states]

        for site in [[1, 0], [0, 1]]:
            for state in states:
                empty_lattice(lattice, state)
                simulation.process_desorb(site, 1)

                if state in states_:
                    self.validate_states(lattice, ('E', 'E', state[2]))

                else:
                    self.validate_states(lattice, state)

        states_ = [(particle, 'O', 'O') for particle in base_states]

        for site in [[1, 2], [2, 1]]:
            for state in states:
                empty_lattice(lattice, state)
                simulation.process_desorb(site, 1)

                if state in states_:
                    self.validate_states(lattice, (state[0], 'E', 'E'))

                else:
                    self.validate_states(lattice, state)

        # Carbon monoxide desorption.
        states_ = [
            ('CO', *particle)
            for particle in product(*([base_states] * 2))
        ]
        site = [0]

        for state in states:
            empty_lattice(lattice, state)
            simulation.process_desorb(site, 4)

            if state in states_:
                self.validate_states(lattice, ('E', *state[1:]))

            else:
                self.validate_states(lattice, state)

        states_ = [
            (particle[0], 'CO', particle[1])
            for particle in product(*([base_states] * 2))
        ]
        site = [1]

        for state in states:
            empty_lattice(lattice, state)
            simulation.process_desorb(site, 4)

            if state in states_:
                self.validate_states(lattice, (state[0], 'E', state[2]))

            else:
                self.validate_states(lattice, state)

        states_ = [
            (*particle, 'CO')
            for particle in product(*([base_states] * 2))
        ]
        site = [2]

        for state in states:
            empty_lattice(lattice, state)
            simulation.process_desorb(site, 4)

            if state in states_:
                self.validate_states(lattice, (*state[:-1], 'E'))

            else:
                self.validate_states(lattice, state)

    def test_process_diffusion(self):
        """
            Tests that the diffusion processes are consistently carried out.
        """
        # ---------------------------------------------------------------------
        # Auxiliary functions.
        # ---------------------------------------------------------------------

        def empty_lattice(lattice0: list, state0: tuple = tuple()) -> None:
            """
                Sets the lattice to the given state. If the state is empty, the
                lattice will be emptied.

                :param lattice0: The lattice whose state is to be set.

                :param state0: The state to which the respective lattice
                 components must be set. Must match the dimensions.
            """
            if state0 == tuple():
                for i, _ in enumerate(lattice0):
                    lattice0[i] = 'E'
                return

            for i, state0_ in enumerate(state0):
                lattice0[i] = state0_

        # ---------------------------------------------------------------------
        # Implementation
        # ---------------------------------------------------------------------

        simulation, _ = TestCOOxidationKMC.get_simulation()
        lattice = simulation.lattice
        base_states = ['CO', 'O', 'E']
        states = list(product(*([base_states] * 3)))

        for particle in [('O', 2), ('CO', 5)]:
            states_ = [
                (particle[0], 'E', particle_)
                for particle_ in base_states
            ]
            site = [0, 1]

            for state in states:
                empty_lattice(lattice, state)
                simulation.process_diffusion(site, particle[1])

                if state in states_:
                    self.validate_states(lattice, ('E', particle[0], state[2]))

                else:
                    self.validate_states(lattice, state)

            states_ = [
                ('E', particle[0], particle_)
                for particle_ in base_states
            ]
            site = [1, 0]

            for state in states:
                empty_lattice(lattice, state)
                simulation.process_diffusion(site, particle[1])

                if state in states_:
                    self.validate_states(lattice, (particle[0], 'E', state[2]))

                else:
                    self.validate_states(lattice, state)

            states_ = [
                (particle_, particle[0], 'E')
                for particle_ in base_states
            ]
            site = [1, 2]

            for state in states:
                empty_lattice(lattice, state)
                simulation.process_diffusion(site, particle[1])

                if state in states_:
                    self.validate_states(lattice, (state[0], 'E', particle[0]))

                else:
                    self.validate_states(lattice, state)

            states_ = [
                (particle_, 'E', particle[0])
                for particle_ in base_states
            ]
            site = [2, 1]

            for state in states:
                empty_lattice(lattice, state)
                simulation.process_diffusion(site, particle[1])

                if state in states_:
                    self.validate_states(lattice, (state[0], particle[0], 'E'))

                else:
                    self.validate_states(lattice, state)

    def test_process_reaction(self):
        """
            Tests that reactions are properly carried out.
        """
        # ---------------------------------------------------------------------
        # Auxiliary functions.
        # ---------------------------------------------------------------------

        def empty_lattice(lattice0: list, state0: tuple = tuple()) -> None:
            """
                Sets the lattice to the given state. If the state is empty, the
                lattice will be emptied.

                :param lattice0: The lattice whose state is to be set.

                :param state0: The state to which the respective lattice
                 components must be set. Must match the dimensions.
            """
            if state0 == tuple():
                for i, _ in enumerate(lattice0):
                    lattice0[i] = 'E'
                return

            for i, state0_ in enumerate(state0):
                lattice0[i] = state0_

        # ---------------------------------------------------------------------
        # Implementation
        # ---------------------------------------------------------------------

        simulation, _ = TestCOOxidationKMC.get_simulation()
        lattice = simulation.lattice
        base_states = ['CO', 'O', 'E']
        states = list(product(*([base_states] * 3)))

        # Langmuir-Hinshelwood.
        states_ = [('O', 'CO', particle) for particle in base_states]
        states_.extend([('CO', 'O', particle) for particle in base_states])

        for site in [[1, 0], [0, 1]]:
            for state in states:
                empty_lattice(lattice, state)
                simulation.process_reaction(site, 6)

                if state in states_:
                    self.validate_states(lattice, ('E', 'E', state[2]))

                else:
                    self.validate_states(lattice, state)

        states_ = [(particle, 'O', 'CO') for particle in base_states]
        states_.extend([(particle, 'CO', 'O') for particle in base_states])

        for site in [[1, 2], [2, 1]]:
            for state in states:
                empty_lattice(lattice, state)
                simulation.process_reaction(site, 6)

                if state in states_:
                    self.validate_states(lattice, (state[0], 'E', 'E'))

                else:
                    self.validate_states(lattice, state)

        # Elay-Rideal.
        states_ = [
            ('O', particle[0], particle[1])
            for particle in product(*([base_states]*2))
        ]
        site = [0]

        for state in states:
            empty_lattice(lattice, state)
            simulation.process_reaction(site, 7)

            if state in states_:
                self.validate_states(lattice, ('E', *state[1:]))

            else:
                self.validate_states(lattice, state)

        states_ = [
            (particle[0], 'O', particle[1])
            for particle in product(*([base_states] * 2))
        ]
        site = [1]

        for state in states:
            empty_lattice(lattice, state)
            simulation.process_reaction(site, 7)

            if state in states_:
                self.validate_states(lattice, (state[0], 'E', state[2]))

            else:
                self.validate_states(lattice, state)

        states_ = [
            (*particle, 'O')
            for particle in product(*([base_states] * 2))
        ]
        site = [2]

        for state in states:
            empty_lattice(lattice, state)
            simulation.process_reaction(site, 7)

            if state in states_:
                self.validate_states(lattice, (*state[:2], 'E'))

            else:
                self.validate_states(lattice, state)

    def test_statistics_record(self):
        """
            Tests that updating the statistics is working.
        """
        simulation, _ = self.get_simulation()

        for i in range(len(simulation)):
            self.assertEqual(simulation.lattice[i], "E")

        simulation.statistics_record()

        for i in range(len(simulation)):
            self.assertEqual(simulation.statistics[i]['E'], 1)

        simulation.lattice[1] = 'CO'
        simulation.statistics_record()

        for i in range(len(simulation)):
            if not i == 1:
                self.assertEqual(simulation.statistics[i]['E'], 2)
                continue

            self.assertEqual(simulation.statistics[i]['E'], 1)
            self.assertEqual(simulation.statistics[i]['CO'], 1)

        simulation.lattice[1] = 'O'
        simulation.statistics_record()

        for i in range(len(simulation)):
            if not i == 1:
                self.assertEqual(simulation.statistics[i]['E'], 3)
                continue

            self.assertEqual(simulation.statistics[i]['E'], 1)
            self.assertEqual(simulation.statistics[i]['CO'], 1)
            self.assertEqual(simulation.statistics[i]['O'], 1)

    def test_update_counters(self):
        """
            Tests that the site validation method is working.
        """
        length = 10
        simulation, _ = self.get_simulation(length=length)

        counter_steps = -1

        for _ in range(30):
            counter_steps = copy.deepcopy(simulation.counter_steps)
            counter_time = copy.deepcopy(simulation.counter_time)
            simulation.update_counters()

            self.assertLess(counter_steps, simulation.counter_steps)
            self.assertLess(counter_time, simulation.counter_time)

        self.assertTrue(counter_steps, 29)

    def test_validate_in_lattice(self):
        """
            Tests that the site validation method is working.
        """
        length = 10
        simulation, _ = self.get_simulation(length=length)

        valid = list(range(len(simulation)))

        for i in valid:
            self.assertTrue(simulation.validate_in_lattice(i))

        not_valid = [i for i in range(len(simulation), 2 * len(simulation))]
        not_valid.extend([i for i in range(-len(simulation), 0)])

        for i in not_valid:
            self.assertFalse(simulation.validate_in_lattice(i))


if __name__ == '__main__':
    # Run the tests.
    unittest.main()
