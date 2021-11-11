""" Contains the tests for the EquationGenerator class."""

# ------------------------------------------------------------------------------
# Imports.
# ------------------------------------------------------------------------------

# Imports: The unittest module.
import unittest

# Imports: Class to be tested.
from coOxidation.Program.Analytic.equation_generator import EquationGenerator

# ------------------------------------------------------------------------------
# Classes.
# ------------------------------------------------------------------------------


class TestEquationGenerator(unittest.TestCase):
    """ Tests that the different functions of the EquationGenerator class for
        the carbon moxide oxidation model.
    """

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Methods Tests.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Get Methods Tests.
    # --------------------------------------------------------------------------

    def test_get_contracted_state(self):
        """ Tests that the _get_contracted_states function is working properly.
        """
        sites_number = 6
        system = EquationGenerator(sites_number)

        # Define the states to be contracted and test the function.
        states = [
            (('CO', 2), ('E', 3), ('CO', 4), ('O', 5)),
            (('CO', 2), ('E', 3), ('O', 4), ('O', 5)),
            (('CO', 2), ('E', 3), ('E', 4), ('O', 5)),
        ]

        # The state that must result from the contraction.
        resulting_state = (('CO', 2), ('E', 3), ('O', 5))
        contracted_state, original_states = system.get_contracted_state(states, 2)
        self.assertEqual(resulting_state, contracted_state)
        self.assertEqual(states, original_states)

        # Not good states to contract.
        states = [
            (('CO', 2), ('E', 3), ('CO', 4), ('O', 5)),
            (('CO', 2), ('E', 3), ('O', 4), ('O', 6)),
            (('CO', 2), ('E', 3), ('E', 4), ('O', 5)),
        ]

        # The state that must result from the contraction.
        resulting_state = tuple()
        contracted_state, original_states = system.get_contracted_state(states, 2)
        self.assertEqual(contracted_state, resulting_state)
        self.assertEqual(states, original_states)

        # Cannot contract the given index.
        states = [
            (('CO', 2), ('E', 3), ('CO', 4), ('O', 5)),
            (('CO', 2), ('E', 3), ('O', 4), ('O', 5)),
            (('CO', 2), ('E', 3), ('E', 4), ('O', 5)),
        ]

        resulting_state = tuple()
        contracted_state, original_states = system.get_contracted_state(states, 1)
        self.assertEqual(contracted_state, resulting_state)
        self.assertEqual(states, original_states)

        # Index contracts to one.
        states = [
            (('CO', 4),), (('O', 4),), (('E', 4),),
        ]

        # The state that must result from the contraction.
        resulting_state = (1,)
        contracted_state, original_states = system.get_contracted_state(states, 0)
        self.assertEqual(contracted_state, resulting_state)
        self.assertEqual(states, original_states)

    def test_get_decay_states(self):
        """ Tests that the _get_decay_states function is working properly.
        """

        sites_number = 6
        system = EquationGenerator(sites_number)

        # Test that the correct dictionary is returned for one-site state.
        state = (('E', 1),)
        decay_dict = {
            'k.O.ads': [],
            'k.O.des': [],
            'k.O.dif': [],
            'k.CO.ads': [(('CO', 1),)],
            'k.CO.des': [],
            'k.CO.dif': [],
            'k.COO.er': [],
            'k.COO.lh': [],
        }

        # Get the dictionary of states.
        states_of_decay = system.get_decay_states(state, system.get_process_functions())
        self.assertEqual(state, states_of_decay[0])
        self.assertEqual(decay_dict, states_of_decay[1])

        # Test that the correct dictionary is returned for two-site states.
        state = (('E', 1), ('O', 2),)
        decay_dict = {
            'k.O.ads': [],
            'k.O.des': [],
            'k.O.dif': [(('O', 1), ('E', 2),)],
            'k.CO.ads': [(('CO', 1), ('O', 2),)],
            'k.CO.des': [],
            'k.CO.dif': [],
            'k.COO.er': [(('E', 1), ('E', 2),)],
            'k.COO.lh': [],
        }

        # Get the dictionary of states.
        states_of_decay = system.get_decay_states(state, system.get_process_functions())
        self.assertEqual(state, states_of_decay[0])
        self.assertEqual(decay_dict, states_of_decay[1])

        # Test that the correct dictionary is returned for three-site states.
        state = (('CO', 1), ('E', 2), ('E', 3))
        decay_dict = {
            'k.O.ads': [(('CO', 1), ('O', 2), ('O', 3))],
            'k.O.des': [],
            'k.O.dif': [],
            'k.CO.ads': [(('CO', 1), ('CO', 2), ('E', 3)), (('CO', 1), ('E', 2), ('CO', 3))],
            'k.CO.des': [(('E', 1), ('E', 2), ('E', 3))],
            'k.CO.dif': [(('E', 1), ('CO', 2), ('E', 3))],
            'k.COO.lh': [],
            'k.COO.er': []
        }

        # Get the dictionary of states.
        states_of_decay = system.get_decay_states(state, system.get_process_functions())
        self.assertEqual(state, states_of_decay[0])
        self.assertEqual(decay_dict, states_of_decay[1])

    def test_get_is_substate(self):
        """ Tests that the _get_is_substate function is working properly.
        """

        sites_number = 6
        system = EquationGenerator(sites_number)

        # Test for consistency.
        state1 = (('CO', 1), ('CO', 2),)
        state2 = (('CO', 1), ('CO', 2), ('E', 3),)

        # State 1 is a substate of state 2.
        self.assertTrue(system.get_is_substate(state1, state2))
        self.assertFalse(system.get_is_substate(state2, state1))

        # Same states with the entries mixed.
        state1 = (('CO', 2), ('CO', 1),)
        state2 = (('E', 3), ('CO', 1), ('CO', 2),)
        self.assertTrue(system.get_is_substate(state1, state2))

        # State 2 is NOT a substate of state 1.
        self.assertFalse(system.get_is_substate(state2, state1))

        # State 1 cannot be a substate of state 2s.
        state1 = (('CO', 1), ('CO', 2),)
        state2 = (('CO', 1), ('CO', 5), ('E', 6),)
        self.assertFalse(system.get_is_substate(state1, state2))
        self.assertFalse(system.get_is_substate(state2, state1))

    def test_get_multiplicity(self):
        """ Tests that the _get_multiplicity function is working properly.
        """

        sites_number = 6
        system = EquationGenerator(sites_number)

        decay_dict0 = {
            'k.O.ads': [(('CO', 1), ('O', 2), ('O', 3)), (('O', 1), ('O', 2), ('O', 3)),
                        (('CO', 1), ('O', 2), ('O', 3))],
            'k.O.des': [],
            'k.O.dif': [],
            'k.CO.ads': [(('CO', 1), ('CO', 2), ('E', 3)), (('CO', 1), ('E', 2), ('CO', 3))],
            'k.CO.des': [(('E', 1), ('E', 2), ('E', 3))],
            'k.CO.dif': [(('E', 1), ('CO', 2)), (('E', 1), ('CO', 2), ('E', 3)), (('E', 1), ('CO', 2), ('E', 3))],
            'k.COO.er': [],
            'k.COO.lh': [],
        }

        decay_dict_resultant0 = {
            'k.O.ads': [((('CO', 1), ('O', 2), ('O', 3)), 2), ((('O', 1), ('O', 2), ('O', 3)), 1)],
            'k.O.des': [],
            'k.O.dif': [],
            'k.CO.ads': [((('CO', 1), ('CO', 2), ('E', 3)), 1), ((('CO', 1), ('E', 2), ('CO', 3)), 1)],
            'k.CO.des': [((('E', 1), ('E', 2), ('E', 3)), 1)],
            'k.CO.dif': [((('E', 1), ('CO', 2)), 1), ((('E', 1), ('CO', 2), ('E', 3)), 2)],
            'k.COO.er': [],
            'k.COO.lh': [],
        }

        decay_dict_resultant1 = system.get_multiplicity(decay_dict0)
        keys = decay_dict_resultant0.keys()
        for key in keys:
            self.assertEqual(len(decay_dict_resultant0[key]), len(decay_dict_resultant1[key]))
            self.assertEqual(set(decay_dict_resultant0[key]), set(decay_dict_resultant1[key]))

    def test_get_numbering(self):
        """ Tests that the _get_numbering function is working properly.
        """

        number_of_sites = 6
        system = EquationGenerator(number_of_sites)

        # Get requested states for a 2 state site.
        state = ('CO', 'CO')
        outcomes_0 = system.get_numbering(state)

        outcomes_1 = [
            (('CO', 1), ('CO', 2),), (('CO', 2), ('CO', 3),),
            (('CO', 3), ('CO', 4),), (('CO', 4), ('CO', 5),),
            (('CO', 5), ('CO', 6),)
        ]

        # Check that all the states are unique.
        self.assertEqual(len(outcomes_0), len(set(outcomes_0)))
        self.assertEqual(len(outcomes_1), len(set(outcomes_1)))
        for outcome_0 in outcomes_0:
            self.assertTrue(outcome_0 in outcomes_1)

        # Get requested states for a 3 state site.
        state = ('CO', 'CO', 'E')
        outcomes_0 = system.get_numbering(state)
        outcomes_1 = [
            (('CO', 1), ('CO', 2), ('E', 3),),
            (('CO', 2), ('CO', 3), ('E', 4),),
            (('CO', 3), ('CO', 4), ('E', 5),),
            (('CO', 4), ('CO', 5), ('E', 6),)
        ]

        # Check that all the states are unique.
        self.assertEqual(len(outcomes_0), len(set(outcomes_0)))
        self.assertEqual(len(outcomes_1), len(set(outcomes_1)))
        for outcome_0 in outcomes_0:
            self.assertTrue(outcome_0 in outcomes_1)

        # Get requested states for a 7 state site; not possible.
        state = ('CO', 'CO', 'E', 'E', 'E', 'E', 'E')
        with self.assertRaises(ValueError):
            system.get_numbering(state)

        state = tuple()
        with self.assertRaises(ValueError):
            system.get_numbering(state)

    def test_get_states(self):
        """ Tests that the _get_states function is working properly.
        """

        number_of_sites = 6
        system = EquationGenerator(number_of_sites)

        # Get requested states for a 2 state site.
        order = 2
        outcomes_0 = system.get_states(order)
        outcomes_1 = [
            ('CO', 'CO',), ('CO', 'O',), ('CO', 'E',),
            ('O', 'CO',), ('O', 'O',), ('O', 'E',),
            ('E', 'CO',), ('E', 'O',), ('E', 'E',)
        ]

        # Check that all the states are unique.
        self.assertEqual(len(outcomes_0), len(set(outcomes_0)))
        self.assertEqual(len(outcomes_1), len(set(outcomes_1)))
        for outcome_0 in outcomes_0:
            self.assertTrue(outcome_0 in outcomes_1)

    def test_get_states_left(self):
        """ Tests that the _get_states_left function is working properly.
        """

        # Define the states and the number of sites.
        number_of_sites = 3
        system = EquationGenerator(number_of_sites)

        # Get lowest level states requested states for zeroth order.
        system.order = 0
        outcomes_0 = system.get_states_left()
        outcomes_1 = [
            (('CO', 1,),), (('CO', 2,),), (('CO', 3,),),
            (('O', 1,),), (('O', 2,),), (('O', 3,),),
            (('E', 1,),), (('E', 2,),), (('E', 3,),)
        ]

        # Check that all the states are unique.
        self.assertEqual(len(outcomes_0), len(set(outcomes_0)))
        self.assertEqual(len(outcomes_1), len(outcomes_0))
        self.assertEqual(len(outcomes_1), len(set(outcomes_1)))

        # Check that all the outcomes are in the list.
        for outcome_0 in outcomes_0:
            self.assertTrue(outcome_0 in outcomes_1)

        # Get lowest level states requested states for first order.
        system.order = 1
        outcomes_0 = system.get_states_left()
        outcomes_1 = [
            (('CO', 1,),), (('CO', 2,),), (('CO', 3,),),
            (('O', 1,),), (('O', 2,),), (('O', 3,),),
            (('E', 1,),), (('E', 2,),), (('E', 3,),)
        ]

        # Check that all the states are unique.
        self.assertEqual(len(outcomes_0), len(set(outcomes_0)))
        self.assertEqual(len(outcomes_1), len(outcomes_0))
        self.assertEqual(len(outcomes_1), len(set(outcomes_1)))

        # Check that all the outcomes are in the list.
        for outcome_0 in outcomes_0:
            self.assertTrue(outcome_0 in outcomes_1)

        # Get lowest level states requested states for second order.
        system.order = 2
        outcomes_0 = system.get_states_left()
        outcomes_1 = [
            (('CO', 1,),), (('CO', 2,),), (('CO', 3,),),
            (('O', 1,),), (('O', 2,),), (('O', 3,),),
            (('E', 1,),), (('E', 2,),), (('E', 3,),),
            (('CO', 1,), ('CO', 2,),), (('CO', 2,), ('CO', 3,),),
            (('CO', 1,), ('O', 2,),), (('CO', 2,), ('O', 3,),),
            (('CO', 1,), ('E', 2,),), (('CO', 2,), ('E', 3,),),
            (('O', 1,), ('CO', 2,),), (('O', 2,), ('CO', 3,),),
            (('O', 1,), ('O', 2,),), (('O', 2,), ('O', 3,),),
            (('O', 1,), ('E', 2,),), (('O', 2,), ('E', 3,),),
            (('E', 1,), ('CO', 2,),), (('E', 2,), ('CO', 3,),),
            (('E', 1,), ('O', 2,),), (('E', 2,), ('O', 3,),),
            (('E', 1,), ('E', 2,),), (('E', 2,), ('E', 3,),)
        ]

        # Check that all the states are unique.
        self.assertEqual(len(outcomes_0), len(set(outcomes_0)))
        self.assertEqual(len(outcomes_1), len(outcomes_0))
        self.assertEqual(len(outcomes_1), len(set(outcomes_1)))

        # Check that all the outcomes are in the list.
        for outcome_0 in outcomes_0:
            self.assertTrue(outcome_0 in outcomes_1)

    def test_get_states_right(self):
        """ Tests that the _get_states_right function is working properly.
        """

        number_of_sites = 3
        system = EquationGenerator(number_of_sites)

        # Get requested states for the first and zeroth order equations.
        system.order = 0
        outcomes_0_0 = system.get_states_right()
        outcomes_0_1 = system.get_states_right()
        outcomes_1 = [
            (('CO', 1,),), (('CO', 2,),), (('CO', 3,),),
            (('O', 1,),), (('O', 2,),), (('O', 3,),),
            (('E', 1,),), (('E', 2,),), (('E', 3,),),
            (('CO', 1,), ('CO', 2,),), (('CO', 2,), ('CO', 3,),),
            (('CO', 1,), ('O', 2,),), (('CO', 2,), ('O', 3,),),
            (('CO', 1,), ('E', 2,),), (('CO', 2,), ('E', 3,),),
            (('O', 1,), ('CO', 2,),), (('O', 2,), ('CO', 3,),),
            (('O', 1,), ('O', 2,),), (('O', 2,), ('O', 3,),),
            (('O', 1,), ('E', 2,),), (('O', 2,), ('E', 3,),),
            (('E', 1,), ('CO', 2,),), (('E', 2,), ('CO', 3,),),
            (('E', 1,), ('O', 2,),), (('E', 2,), ('O', 3,),),
            (('E', 1,), ('E', 2,),), (('E', 2,), ('E', 3,),)
        ]

        # Check that all the states are unique.
        self.assertEqual(len(outcomes_0_0), len(set(outcomes_0_0)))
        self.assertEqual(len(outcomes_0_1), len(set(outcomes_0_1)))
        self.assertEqual(len(outcomes_1), len(outcomes_0_0))
        self.assertEqual(len(outcomes_1), len(outcomes_0_1))
        self.assertEqual(len(outcomes_1), len(set(outcomes_1)))

        # Check that all the outcomes are in the list.
        for j, outcome_0 in enumerate(outcomes_0_0):
            self.assertTrue(outcome_0 in outcomes_1)
            self.assertTrue(outcomes_0_1[j] in outcomes_1)

        # Get requested states for the exact equations.
        system.order = 34
        outcomes_0_0 = system.get_states_right()

        system.order = 3
        outcomes_0_1 = system.get_states_right()

        # These are the states that must come out.
        outcomes_1 = [
            (('CO', 1,), ('CO', 2,), ('CO', 3,),),
            (('CO', 1,), ('CO', 2,), ('O', 3,),),
            (('CO', 1,), ('CO', 2,), ('E', 3,),),
            (('CO', 1,), ('O', 2,), ('CO', 3,),),
            (('CO', 1,), ('O', 2,), ('O', 3,),),
            (('CO', 1,), ('O', 2,), ('E', 3,),),
            (('CO', 1,), ('E', 2,), ('CO', 3,),),
            (('CO', 1,), ('E', 2,), ('O', 3,),),
            (('CO', 1,), ('E', 2,), ('E', 3,),),
            (('O', 1,), ('CO', 2,), ('CO', 3,),),
            (('O', 1,), ('CO', 2,), ('O', 3,),),
            (('O', 1,), ('CO', 2,), ('E', 3,),),
            (('O', 1,), ('O', 2,), ('CO', 3,),),
            (('O', 1,), ('O', 2,), ('O', 3,),),
            (('O', 1,), ('O', 2,), ('E', 3,),),
            (('O', 1,), ('E', 2,), ('CO', 3,),),
            (('O', 1,), ('E', 2,), ('O', 3,),),
            (('O', 1,), ('E', 2,), ('E', 3,),),
            (('E', 1,), ('CO', 2,), ('CO', 3,),),
            (('E', 1,), ('CO', 2,), ('O', 3,),),
            (('E', 1,), ('CO', 2,), ('E', 3,),),
            (('E', 1,), ('O', 2,), ('CO', 3,),),
            (('E', 1,), ('O', 2,), ('O', 3,),),
            (('E', 1,), ('O', 2,), ('E', 3,),),
            (('E', 1,), ('E', 2,), ('CO', 3,),),
            (('E', 1,), ('E', 2,), ('O', 3,),),
            (('E', 1,), ('E', 2,), ('E', 3,),)
        ]

        # Check that all the states are unique.
        self.assertEqual(len(outcomes_0_0), len(set(outcomes_0_0)))
        self.assertEqual(len(outcomes_0_1), len(set(outcomes_0_1)))
        self.assertEqual(len(outcomes_1), len(outcomes_0_0))
        self.assertEqual(len(outcomes_1), len(outcomes_0_1))
        self.assertEqual(len(outcomes_1), len(set(outcomes_1)))

        # Check that all the outcomes are in the list.
        for j, outcome_0 in enumerate(outcomes_0_0):
            self.assertTrue(outcome_0 in outcomes_1)
            self.assertTrue(outcomes_0_1[j] in outcomes_1)

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # General Properties Tests.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Model Consistency Tests.
    # --------------------------------------------------------------------------

    def test_consistency(self):
        """ Tests that the model is consistent, i.e., the number of rate
            constants must be the same as the number of rates for the processes,
            that must be the same as the number of functions to be applied.
        """

        number_of_sites = 3
        system = EquationGenerator(number_of_sites)

        # Check the rate constants.
        processes_rates = system.get_process_rates()
        self.assertEqual(len(processes_rates), len(set(processes_rates)))

        # Check the process orders.
        process_dict = system.get_process_orders()
        process_orders = set(process_dict[process_rate] for process_rate in processes_rates)
        self.assertTrue(all(map(lambda x: x > 0, process_orders)))
        self.assertGreaterEqual(system.sites_number, max(process_orders))

        # Check the functions are unique.
        process_dict = system.get_process_functions()
        process_functions = set(process_dict[process_rate] for process_rate in processes_rates)
        self.assertEqual(len(process_functions), len(set(process_functions)))

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Processes Tests
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Carbon Monoxide Processes Tests.
    # --------------------------------------------------------------------------

    def test_carbon_monoxide_adsorb(self):
        """ Tests that the _carbon_monoxide_adsorb function is working properly.
        """

        number_of_sites = 3
        system = EquationGenerator(number_of_sites)

        # Adsorption of carbon monoxide, with CO particles.
        mock_state = (('CO', 1), ('E', 2), ('E', 3),)
        outcomes_0 = (('CO', 1), ('CO', 2), ('E', 3),), (('CO', 1), ('E', 2), ('CO', 3),)
        outcomes_1 = system.process_carbon_monoxide_adsorb(mock_state)
        self.assertEqual(len(outcomes_0), len(outcomes_1))
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # Adsorption of carbon monoxide, with O particles.
        mock_state = (('O', 1), ('E', 2), ('E', 3),)
        outcomes_0 = (('O', 1), ('CO', 2), ('E', 3),), (('O', 1), ('E', 2), ('CO', 3),)
        outcomes_1 = system.process_carbon_monoxide_adsorb(mock_state)
        self.assertEqual(len(outcomes_0), len(outcomes_1))
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # NO adsorption of carbon monoxide.
        mock_state = (('O', 1), ('CO', 2), ('O', 3),)
        outcomes_1 = system.process_carbon_monoxide_adsorb(mock_state)
        self.assertEqual(len(outcomes_1), 0)

    def test_carbon_monoxide_desorb(self):
        """ Tests that the _carbon_monoxide_desorb function is working properly.
        """

        number_of_sites = 3
        system = EquationGenerator(number_of_sites)

        # Desorption of carbon monoxide, with CO particles.
        mock_state = (('CO', 1), ('E', 2), ('CO', 3),)
        outcomes_0 = (('CO', 1), ('E', 2), ('E', 3),), (('E', 1), ('E', 2), ('CO', 3),)
        outcomes_1 = system.process_carbon_monoxide_desorb(mock_state)
        self.assertEqual(len(outcomes_0), len(outcomes_1))
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # Desorption of carbon monoxide, with O particles.
        mock_state = (('O', 1), ('E', 2), ('CO', 3),)
        outcomes_0 = (('O', 1), ('E', 2), ('E', 3),),
        outcomes_1 = system.process_carbon_monoxide_desorb(mock_state)
        self.assertEqual(len(outcomes_0), len(outcomes_1))
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # NO desorption of carbon monoxide.
        mock_state = (('O', 1), ('E', 2), ('O', 3),)
        outcomes_1 = system.process_carbon_monoxide_desorb(mock_state)
        self.assertEqual(len(outcomes_1), 0)

    def test_carbon_monoxide_diffusion(self):
        """ Tests that the _carbon_monoxide_diffusion function is working properly.
        """

        number_of_sites = 6
        system = EquationGenerator(number_of_sites)

        # Diffusion of carbon monoxide, with CO particles.
        mock_state = (('CO', 1), ('E', 2), ('CO', 3), ('E', 6),)
        outcomes_0 = (('E', 1), ('CO', 2), ('CO', 3), ('E', 6),), (('CO', 1), ('CO', 2), ('E', 3), ('E', 6),)
        outcomes_1 = system.process_carbon_monoxide_diffusion(mock_state)
        self.assertEqual(len(outcomes_0), len(outcomes_1))
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # Diffusion of carbon monoxide, with O particles.
        mock_state = (('O', 1), ('E', 2), ('CO', 3),)
        outcomes_0 = (('O', 1), ('CO', 2), ('E', 3),),
        outcomes_1 = system.process_carbon_monoxide_diffusion(mock_state)
        self.assertEqual(len(outcomes_0), len(outcomes_1))
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # NO adsorption of carbon monoxide.
        mock_state = (('O', 1), ('E', 2), ('O', 3),)
        outcomes_1 = system.process_carbon_monoxide_diffusion(mock_state)
        self.assertEqual(len(outcomes_1), 0)

        # Diffusion of carbon monoxide with only one state.
        mock_state = (('E', 1),)
        outcomes_1 = system.process_carbon_monoxide_diffusion(mock_state)
        self.assertEqual(len(outcomes_1), 0)

    # --------------------------------------------------------------------------
    # Oxygen Processes Tests.
    # --------------------------------------------------------------------------

    def test_oxygen_adsorb(self):
        """ Tests that the _oxygen_adsorb function is working properly.
        """

        number_of_sites = 6
        system = EquationGenerator(number_of_sites)

        # Adsorption of oxygen, with CO particles.
        mock_state = (('CO', 1), ('E', 2), ('E', 3), ('E', 6))
        outcomes_0 = (('CO', 1), ('O', 2), ('O', 3), ('E', 6)),
        outcomes_1 = system.process_oxygen_adsorb(mock_state)
        self.assertEqual(len(outcomes_0), len(outcomes_1))
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # Adsorption of oxygen, with O particles.
        mock_state = (('E', 1), ('E', 2), ('O', 3), ('E', 4),)
        outcomes_0 = (('O', 1), ('O', 2), ('O', 3), ('E', 4),),
        outcomes_1 = system.process_oxygen_adsorb(mock_state)
        self.assertEqual(len(outcomes_0), len(outcomes_1))
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # NO adsorption of oxygen.
        mock_state = (('O', 1), ('E', 2), ('O', 3), ('E', 4),)
        outcomes_1 = system.process_oxygen_adsorb(mock_state)
        self.assertEqual(len(outcomes_1), 0)

    def test_oxygen_desorb(self):
        """ Tests that the _oxygen_desorb function is working properly.
        """

        number_of_sites = 6
        system = EquationGenerator(number_of_sites)

        # Desorption of oxygen, with CO particles.
        mock_state = (('CO', 1), ('O', 2), ('O', 3), ('O', 6))
        outcomes_0 = (('CO', 1), ('E', 2), ('E', 3), ('O', 6)),
        outcomes_1 = system.process_oxygen_desorb(mock_state)
        self.assertEqual(len(outcomes_0), len(outcomes_1))
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # Desorption of oxygen, with O particles.
        mock_state = (('E', 1), ('O', 2), ('O', 3), ('O', 4),)
        outcomes_0 = (('E', 1), ('E', 2), ('E', 3), ('O', 4),), (('E', 1), ('O', 2), ('E', 3), ('E', 4),)
        outcomes_1 = system.process_oxygen_desorb(mock_state)
        self.assertEqual(len(outcomes_0), len(outcomes_1))
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # NO desorption of oxygen.
        mock_state = (('O', 1), ('E', 2), ('O', 3), ('O', 5),)
        outcomes_1 = system.process_oxygen_desorb(mock_state)
        self.assertEqual(len(outcomes_1), 0)

    def test_oxygen_diffusion(self):
        """ Tests that the _oxygen_diffusion function is working properly.
        """

        number_of_sites = 6
        system = EquationGenerator(number_of_sites)

        # Desorption of oxygen, with CO particles.
        mock_state = (('CO', 1), ('O', 2), ('E', 3), ('O', 6))
        outcomes_0 = (('CO', 1), ('E', 2), ('O', 3), ('O', 6)),
        outcomes_1 = system.process_oxygen_diffusion(mock_state)
        self.assertEqual(len(outcomes_0), len(outcomes_1))
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # Desorption of oxygen, with O particles.
        mock_state = (('E', 1), ('O', 2), ('E', 3), ('O', 4),)
        outcomes_0 = (
            (('O', 1), ('E', 2), ('E', 3), ('O', 4),),
            (('E', 1), ('E', 2), ('O', 3), ('O', 4),),
            (('E', 1), ('O', 2), ('O', 3), ('E', 4),)
        )
        outcomes_1 = system.process_oxygen_diffusion(mock_state)
        self.assertEqual(len(outcomes_0), len(outcomes_1))
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # NO desorption of oxygen.
        mock_state = (('O', 1), ('CO', 2), ('O', 3), ('O', 5),)
        outcomes_1 = system.process_oxygen_diffusion(mock_state)
        self.assertEqual(len(outcomes_1), 0)


if __name__ == '__main__':
    # Run the tests.
    unittest.main()
