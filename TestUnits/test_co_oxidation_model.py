""" Contains the tests for the COOxidationEquationGenerator class.
"""
import os
path1 = os.getcwd()
os.chdir("..")

# Imports: The unittest module.
import unittest

# Imports: Class to be tested.
from Program.co_oxidation_model import COOxidationEquationGenerator as Co


class TestCOOxidationEquationGenerator(unittest.TestCase):
    """ Tests that the different functions of the COOxidationEquationGenerator
        class are working properly.
    """

    # --------------------------------------------------------------------------
    # Model Consistency Tests.
    # --------------------------------------------------------------------------

    def test_consistency(self):
        """ Tests that the model is consistent, i.e., the number of rate
            constants must be the same as the number of rates for the processes,
            that must be the same as the number of functions to be applied.
        """

        # Define the states and the number of sites.
        states = {'E', 'CO', 'O'}
        number_of_sites = 3

        # Create a CO oxidation system with three sites.
        system = Co(number_of_sites, states)

        # ----------------------------------------------------------------------
        # Check the rate constants.
        # ----------------------------------------------------------------------

        # Get the rate constants.
        processes_rates = system._get_process_rates()

        # Make sure they are unique.
        self.assertEqual(len(processes_rates), len(set(processes_rates)))

        # ----------------------------------------------------------------------
        # Check the process orders.
        # ----------------------------------------------------------------------

        # Get the UNIQUE process orders.
        process_dict = system._get_process_orders()
        process_orders = set(process_dict[process_rate] for process_rate in processes_rates)

        # Make sure they are all positive.
        self.assertTrue(all(map(lambda x: x > 0, process_orders)))

        # Make sure the maximum order is less than the number of sites.
        self.assertGreaterEqual(system.number_of_sites, max(process_orders))

        # ----------------------------------------------------------------------
        # Check the functions are unique.
        # ----------------------------------------------------------------------

        # Get the UNIQUE process orders.
        process_dict = system._get_process_functions()
        process_functions = set(process_dict[process_rate] for process_rate in processes_rates)

        # Make sure they are all unique.
        self.assertEqual(len(process_functions), len(set(process_functions)))

    # --------------------------------------------------------------------------
    # Carbon Monoxide Processes Tests.
    # --------------------------------------------------------------------------

    def test_carbon_monoxide_adsorb(self):
        """ Tests that the _carbon_monoxide_adsorb function is working properly.
        """

        # Define the states and the number of sites.
        states = {'E', 'CO', 'O'}
        number_of_sites = 3

        # Create a CO oxidation system with three sites.
        system = Co(number_of_sites, states)

        # ----------------------------------------------------------------------
        # Adsorption of carbon monoxide, with CO particles.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('CO', 1), ('E', 2), ('E', 3),)

        # The list with the possible outcome states.
        outcomes_0 = (('CO', 1), ('CO', 2), ('E', 3),), (('CO', 1), ('E', 2), ('CO', 3),)

        # Apply the carbon monoxide adsorption operation to the state.
        outcomes_1 = system._carbon_monoxide_adsorb(mock_state)

        # Check that the outcomes have the same length.
        self.assertEqual(len(outcomes_0), len(outcomes_1))

        # Check that all the requested states are in the outcome state.
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # ----------------------------------------------------------------------
        # Adsorption of carbon monoxide, with O particles.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('O', 1), ('E', 2), ('E', 3),)

        # The list with the possible outcome states.
        outcomes_0 = (('O', 1), ('CO', 2), ('E', 3),), (('O', 1), ('E', 2), ('CO', 3),)

        # Apply the carbon monoxide adsorption operation to the state.
        outcomes_1 = system._carbon_monoxide_adsorb(mock_state)

        # Check that the outcomes have the same length.
        self.assertEqual(len(outcomes_0), len(outcomes_1))

        # Check that all the requested states are in the outcome state.
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # ----------------------------------------------------------------------
        # NO adsorption of carbon monoxide.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('O', 1), ('CO', 2), ('O', 3),)

        # Apply the carbon monoxide adsorption operation to the state.
        outcomes_1 = system._carbon_monoxide_adsorb(mock_state)

        # Check that there NO outcomes.
        self.assertEqual(len(outcomes_1), 0)

    def test_carbon_monoxide_desorb(self):
        """ Tests that the _carbon_monoxide_desorb function is working properly.
        """

        # Define the states and the number of sites.
        states = {'E', 'CO', 'O'}
        number_of_sites = 3

        # Create a CO oxidation system with three sites.
        system = Co(number_of_sites, states)

        # ----------------------------------------------------------------------
        # Desorption of carbon monoxide, with CO particles.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('CO', 1), ('E', 2), ('CO', 3),)

        # The list with the possible outcome states.
        outcomes_0 = (('CO', 1), ('E', 2), ('E', 3),), (('E', 1), ('E', 2), ('CO', 3),)

        # Apply the carbon monoxide desorption operation to the state.
        outcomes_1 = system._carbon_monoxide_desorb(mock_state)

        # Check that the outcomes have the same length.
        self.assertEqual(len(outcomes_0), len(outcomes_1))

        # Check that all the requested states are in the outcome state.
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # ----------------------------------------------------------------------
        # Desorption of carbon monoxide, with O particles.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('O', 1), ('E', 2), ('CO', 3),)

        # The list with the possible outcome states.
        outcomes_0 = (('O', 1), ('E', 2), ('E', 3),),

        # Apply the carbon monoxide desorption operation to the state.
        outcomes_1 = system._carbon_monoxide_desorb(mock_state)

        # Check that the outcomes have the same length.
        self.assertEqual(len(outcomes_0), len(outcomes_1))

        # Check that all the requested states are in the outcome state.
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # ----------------------------------------------------------------------
        # NO desorption of carbon monoxide.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('O', 1), ('E', 2), ('O', 3),)

        # Apply the carbon monoxide desorption operation to the state.
        outcomes_1 = system._carbon_monoxide_desorb(mock_state)

        # Check that there NO outcomes.
        self.assertEqual(len(outcomes_1), 0)

    def test_carbon_monoxide_diffusion(self):
        """ Tests that the _carbon_monoxide_diffusion function is working properly.
        """

        # Define the states and the number of sites.
        states = {'E', 'CO', 'O'}
        number_of_sites = 6

        # Create a CO oxidation system with six sites.
        system = Co(number_of_sites, states)

        # ----------------------------------------------------------------------
        # Diffusion of carbon monoxide, with CO particles.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('CO', 1), ('E', 2), ('CO', 3), ('E', 6),)

        # The list with the possible outcome states.
        outcomes_0 = (('E', 1), ('CO', 2), ('CO', 3), ('E', 6),), (('CO', 1), ('CO', 2), ('E', 3), ('E', 6),)

        # Apply the carbon monoxide diffusion operation to the state.
        outcomes_1 = system._carbon_monoxide_diffusion(mock_state)

        # Check that the outcomes have the same length.
        self.assertEqual(len(outcomes_0), len(outcomes_1))

        # Check that all the requested states are in the outcome state.
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # ----------------------------------------------------------------------
        # diffusion of carbon monoxide, with O particles.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('O', 1), ('E', 2), ('CO', 3),)

        # The list with the possible outcome states.
        outcomes_0 = (('O', 1), ('CO', 2), ('E', 3),),

        # Apply the carbon monoxide diffusion operation to the state.
        outcomes_1 = system._carbon_monoxide_diffusion(mock_state)

        # Check that the outcomes have the same length.
        self.assertEqual(len(outcomes_0), len(outcomes_1))

        # Check that all the requested states are in the outcome state.
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # ----------------------------------------------------------------------
        # NO adsorption of carbon monoxide.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('O', 1), ('E', 2), ('O', 3),)

        # Apply the carbon monoxide diffusion operation to the state.
        outcomes_1 = system._carbon_monoxide_diffusion(mock_state)

        # Check that there NO outcomes.
        self.assertEqual(len(outcomes_1), 0)

    # --------------------------------------------------------------------------
    # Oxygen Processes Tests.
    # --------------------------------------------------------------------------

    def test_oxygen_adsorb(self):
        """ Tests that the _oxygen_adsorb function is working properly.
        """

        # Define the states and the number of sites.
        states = {'E', 'CO', 'O'}
        number_of_sites = 6

        # Create a CO oxidation system with six sites.
        system = Co(number_of_sites, states)

        # ----------------------------------------------------------------------
        # Adsorption of oxygen, with CO particles.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('CO', 1), ('E', 2), ('E', 3), ('E', 6))

        # The list with the possible outcome states.
        outcomes_0 = (('CO', 1), ('O', 2), ('O', 3), ('E', 6)),

        # Apply the oxygen adsorption operation to the state.
        outcomes_1 = system._oxygen_adsorb(mock_state)

        # Check that the outcomes have the same length.
        self.assertEqual(len(outcomes_0), len(outcomes_1))

        # Check that all the requested states are in the outcome state.
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # ----------------------------------------------------------------------
        # Adsorption of oxygen, with O particles.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('E', 1), ('E', 2), ('O', 3), ('E', 4),)

        # The list with the possible outcome states.
        outcomes_0 = (('O', 1), ('O', 2), ('O', 3), ('E', 4),),

        # Apply the oxygen adsorption operation to the state.
        outcomes_1 = system._oxygen_adsorb(mock_state)

        # Check that the outcomes have the same length.
        self.assertEqual(len(outcomes_0), len(outcomes_1))

        # Check that all the requested states are in the outcome state.
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # ----------------------------------------------------------------------
        # NO adsorption of oxygen.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('O', 1), ('E', 2), ('O', 3), ('E', 4),)

        # Apply the oxygen adsorption operation to the state.
        outcomes_1 = system._oxygen_adsorb(mock_state)

        # Check that there NO outcomes.
        self.assertEqual(len(outcomes_1), 0)

    def test_oxygen_desorb(self):
        """ Tests that the _oxygen_desorb function is working properly.
        """

        # Define the states and the number of sites.
        states = {'E', 'CO', 'O'}
        number_of_sites = 6

        # Create a CO oxidation system with six sites.
        system = Co(number_of_sites, states)

        # ----------------------------------------------------------------------
        # Desorption of oxygen, with CO particles.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('CO', 1), ('O', 2), ('O', 3), ('O', 6))

        # The list with the possible outcome states.
        outcomes_0 = (('CO', 1), ('E', 2), ('E', 3), ('O', 6)),

        # Apply the oxygen desorption operation to the state.
        outcomes_1 = system._oxygen_desorb(mock_state)

        # Check that the outcomes have the same length.
        self.assertEqual(len(outcomes_0), len(outcomes_1))

        # Check that all the requested states are in the outcome state.
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # ----------------------------------------------------------------------
        # Desorption of oxygen, with O particles.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('E', 1), ('O', 2), ('O', 3), ('O', 4),)

        # The list with the possible outcome states.
        outcomes_0 = (('E', 1), ('E', 2), ('E', 3), ('O', 4),), (('E', 1), ('O', 2), ('E', 3), ('E', 4),)

        # Apply the oxygen desorption operation to the state.
        outcomes_1 = system._oxygen_desorb(mock_state)

        # Check that the outcomes have the same length.
        self.assertEqual(len(outcomes_0), len(outcomes_1))

        # Check that all the requested states are in the outcome state.
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # ----------------------------------------------------------------------
        # NO desorption of oxygen.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('O', 1), ('E', 2), ('O', 3), ('O', 5),)

        # Apply the oxygen desorption operation to the state.
        outcomes_1 = system._oxygen_desorb(mock_state)

        # Check that there NO outcomes.
        self.assertEqual(len(outcomes_1), 0)

    def test_oxygen_diffusion(self):
        """ Tests that the _oxygen_diffusion function is working properly.
        """

        # Define the states and the number of sites.
        states = {'E', 'CO', 'O'}
        number_of_sites = 6

        # Create a CO oxidation system with six sites.
        system = Co(number_of_sites, states)

        # ----------------------------------------------------------------------
        # Desorption of oxygen, with CO particles.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('CO', 1), ('O', 2), ('E', 3), ('O', 6))

        # The list with the possible outcome states.
        outcomes_0 = (('CO', 1), ('E', 2), ('O', 3), ('O', 6)),

        # Apply the oxygen diffusion operation to the state.
        outcomes_1 = system._oxygen_diffusion(mock_state)

        # Check that the outcomes have the same length.
        self.assertEqual(len(outcomes_0), len(outcomes_1))

        # Check that all the requested states are in the outcome state.
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # ----------------------------------------------------------------------
        # Desorption of oxygen, with O particles.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('E', 1), ('O', 2), ('E', 3), ('O', 4),)

        # The list with the possible outcome states.
        outcomes_0 = (('O', 1), ('E', 2), ('E', 3), ('O', 4),), (('E', 1), ('E', 2), ('O', 3), ('O', 4),),\
                     (('E', 1), ('O', 2), ('O', 3), ('E', 4),)

        # Apply the oxygen diffusion operation to the state.
        outcomes_1 = system._oxygen_diffusion(mock_state)

        # Check that the outcomes have the same length.
        self.assertEqual(len(outcomes_0), len(outcomes_1))

        # Check that all the requested states are in the outcome state.
        for i, outcome_1 in enumerate(outcomes_1):
            self.assertTrue(outcome_1 in outcomes_0)

        # ----------------------------------------------------------------------
        # NO desorption of oxygen.
        # ----------------------------------------------------------------------

        # Define a valid mock state.
        mock_state = (('O', 1), ('CO', 2), ('O', 3), ('O', 5),)

        # Apply the oxygen diffusion operation to the state.
        outcomes_1 = system._oxygen_diffusion(mock_state)

        # Check that there NO outcomes.
        self.assertEqual(len(outcomes_1), 0)

    # --------------------------------------------------------------------------
    # Other Methods Tests.
    # --------------------------------------------------------------------------

    def test_get_contracted_state(self):
        """ Tests that the _get_contracted_states function is working properly.
        """

        # Define the states and the number of sites.
        states = {'E', 'CO', 'O'}
        number_of_sites = 6

        # Create a CO oxidation system with six sites.
        system = Co(number_of_sites, states)

        # ----------------------------------------------------------------------
        # Define the states to be contracted and test the function.
        # ----------------------------------------------------------------------

        # Define the collection of states to be contracted.
        states = (
            (('CO', 2), ('E', 3), ('CO', 4), ('O', 5)),
            (('CO', 2), ('E', 3), ('O', 4), ('O', 5)),
            (('CO', 2), ('E', 3), ('E', 4), ('O', 5)),
        )

        # The state that must result from the contraction.
        resulting_state = (('CO', 2), ('E', 3), ('O', 5))

        # Apply the contraction function to the second index.
        contracted_state, original_states = system._get_contracted_state(states, 2)

        # The resulting state must be the same as the contracted state.
        self.assertTupleEqual(resulting_state, contracted_state)
        self.assertTupleEqual(states, original_states)

        # ----------------------------------------------------------------------
        # Not good states to contract.
        # ----------------------------------------------------------------------

        # Define the collection of states to be contracted.
        states = (
            (('CO', 2), ('E', 3), ('CO', 4), ('O', 5)),
            (('CO', 2), ('E', 3), ('O', 4), ('O', 6)),
            (('CO', 2), ('E', 3), ('E', 4), ('O', 5)),
        )

        # The state that must result from the contraction.
        resulting_state = tuple([])

        # Apply the contraction function to the second index.
        contracted_state, original_states = system._get_contracted_state(states, 2)

        # The resulting state must be the same as the contracted state.
        self.assertTupleEqual(contracted_state, resulting_state)
        self.assertTupleEqual(states, original_states)

        # ----------------------------------------------------------------------
        # Cannot contract the given index.
        # ----------------------------------------------------------------------

        # Define the collection of states to be contracted.
        states = (
            (('CO', 2), ('E', 3), ('CO', 4), ('O', 5)),
            (('CO', 2), ('E', 3), ('O', 4), ('O', 5)),
            (('CO', 2), ('E', 3), ('E', 4), ('O', 5)),
        )

        # The state that must result from the contraction.
        resulting_state = tuple([])

        # Apply the contraction function to the second index.
        contracted_state, original_states = system._get_contracted_state(states, 1)

        # The resulting state must be the same as the contracted state.
        self.assertTupleEqual(contracted_state, resulting_state)
        self.assertTupleEqual(states, original_states)

        # ----------------------------------------------------------------------
        # Index contracts to one.
        # ----------------------------------------------------------------------

        # Define the collection of states to be contracted.
        states = (
            (('CO', 4),),
            (('O', 4),),
            (('E', 4),),
        )

        # The state that must result from the contraction.
        resulting_state = (1,)

        # Apply the contraction function to the second index.
        contracted_state, original_states = system._get_contracted_state(states, 0)

        # The resulting state must be the same as the contracted state.
        self.assertTupleEqual(contracted_state, resulting_state)
        self.assertTupleEqual(states, original_states)

    def test_get_decay_states(self):
        """ Tests that the _get_decay_states function is working properly.
        """

        # Define the states and the number of sites.
        states = {'E', 'CO', 'O'}
        number_of_sites = 6

        # Create a CO oxidation system with six sites.
        system = Co(number_of_sites, states)

        # ----------------------------------------------------------------------
        # Test that the correct dictionary is returned for one-site state.
        # ----------------------------------------------------------------------

        # Define the state that will be operated on by the different processes.
        state = (('E', 1),)

        # Set the decay dictionary.
        decay_dict = {
            'k.O.ads':  [],
            'k.O.des':  [],
            'k.O.dif':  [],
            'k.CO.ads': [(('CO', 1),)],
            'k.CO.des': [],
            'k.CO.dif': [],
            'k.COO.lh': [],
            'k.COO.el': []
        }

        # Get the dictionary of states.
        states_of_decay = system._get_decay_states(state, system._get_process_functions())

        # Original state must also be returned.
        self.assertEqual(state, states_of_decay[0])

        # Dictionaries must be the same.
        self.assertEqual(decay_dict, states_of_decay[1])

        # ----------------------------------------------------------------------
        # Test that the correct dictionary is returned for two-site states.
        # ----------------------------------------------------------------------

        # Define the state that will be operated on by the different processes.
        state = (('E', 1), ('O', 2),)

        # Set the decay dictionary.
        decay_dict = {
            'k.O.ads': [],
            'k.O.des': [],
            'k.O.dif': [(('O', 1), ('E', 2),)],
            'k.CO.ads': [(('CO', 1), ('O', 2),)],
            'k.CO.des': [],
            'k.CO.dif': [],
            'k.COO.lh': [],
            'k.COO.el': [(('E', 1), ('E', 2),)]
        }

        # Get the dictionary of states.
        states_of_decay = system._get_decay_states(state, system._get_process_functions())

        # Original state must also be returned.
        self.assertEqual(state, states_of_decay[0])

        # Dictionaries must be the same.
        self.assertEqual(decay_dict, states_of_decay[1])

        # ----------------------------------------------------------------------
        # Test that the correct dictionary is returned for three-site states.
        # ----------------------------------------------------------------------

        # Define the state that will be operated on by the different processes.
        state = (('CO', 1), ('E', 2), ('E', 3))

        # Set the decay dictionary.
        decay_dict = {
            'k.O.ads':  [(('CO', 1), ('O', 2), ('O', 3))],
            'k.O.des':  [],
            'k.O.dif':  [],
            'k.CO.ads': [(('CO', 1), ('CO', 2), ('E', 3)), (('CO', 1), ('E', 2), ('CO', 3))],
            'k.CO.des': [(('E', 1), ('E', 2), ('E', 3))],
            'k.CO.dif': [(('E', 1), ('CO', 2), ('E', 3))],
            'k.COO.lh': [],
            'k.COO.el': []
        }

        # Get the dictionary of states.
        states_of_decay = system._get_decay_states(state, system._get_process_functions())

        # Original state must also be returned.
        self.assertEqual(state, states_of_decay[0])

        # Dictionaries must be the same.
        self.assertEqual(decay_dict, states_of_decay[1])

    def test_get_is_substate(self):
        """ Tests that the _get_is_substate function is working properly.
        """

        # Define the states and the number of sites.
        states = {'E', 'CO', 'O'}
        number_of_sites = 6

        # Create a CO oxidation system with six sites.
        system = Co(number_of_sites, states)

        # ----------------------------------------------------------------------
        # Test for consistency.
        # ----------------------------------------------------------------------

        # Get the states to be examined.
        state1 = (('CO', 1), ('CO', 2),)
        state2 = (('CO', 1), ('CO', 2), ('E', 3),)

        # State 1 is a substate of state 2.
        self.assertTrue(system._get_is_substate(state1, state2))

        # State 2 is NOT a substate of state 1.
        self.assertFalse(system._get_is_substate(state2, state1))

        # ----------------------------------------------------------------------
        # Same states with the entries mixed.
        # ----------------------------------------------------------------------

        # Get the states to be examined.
        state1 = (('CO', 2), ('CO', 1),)
        state2 = (('E', 3), ('CO', 1), ('CO', 2),)

        # State 1 is a substate of state 2.
        self.assertTrue(system._get_is_substate(state1, state2))

        # State 2 is NOT a substate of state 1.
        self.assertFalse(system._get_is_substate(state2, state1))

        # ----------------------------------------------------------------------
        # State 1 cannot be a substate of state 2s.
        # ----------------------------------------------------------------------

        # Get the states to be examined.
        state1 = (('CO', 1), ('CO', 2), )
        state2 = (('CO', 1), ('CO', 5), ('E', 6),)

        # State 1 is NOT substate of state 2.
        self.assertFalse(system._get_is_substate(state1, state2))

        # State 2 is NOT a substate of state 1.
        self.assertFalse(system._get_is_substate(state2, state1))

    def test_get_multiplicity(self):
        """ Tests that the _get_multiplicity function is working properly.
        """

        # Define the states and the number of sites.
        states = {'E', 'CO', 'O'}
        number_of_sites = 6

        # Create a CO oxidation system with six sites.
        system = Co(number_of_sites, states)

        # ----------------------------------------------------------------------
        # Test that the correct dictionary is returned for three-site states.
        # ----------------------------------------------------------------------

        # Set the (fictitious) decay dictionary.
        decay_dict0 = {
            'k.O.ads': [(('CO', 1), ('O', 2), ('O', 3)), (('O', 1), ('O', 2), ('O', 3)), (('CO', 1), ('O', 2), ('O', 3))],
            'k.O.des': [],
            'k.O.dif': [],
            'k.CO.ads': [(('CO', 1), ('CO', 2), ('E', 3)), (('CO', 1), ('E', 2), ('CO', 3))],
            'k.CO.des': [(('E', 1), ('E', 2), ('E', 3))],
            'k.CO.dif': [(('E', 1), ('CO', 2)), (('E', 1), ('CO', 2), ('E', 3)), (('E', 1), ('CO', 2), ('E', 3))],
            'k.COO.lh': [],
            'k.COO.el': []
        }

        decay_dict_resultant0 = {
            'k.O.ads': [((('CO', 1), ('O', 2), ('O', 3)), 2), ((('O', 1), ('O', 2), ('O', 3)),1)],
            'k.O.des': [],
            'k.O.dif': [],
            'k.CO.ads': [((('CO', 1), ('CO', 2), ('E', 3)) , 1), ((('CO', 1), ('E', 2), ('CO', 3)), 1)],
            'k.CO.des': [((('E', 1), ('E', 2), ('E', 3)), 1)],
            'k.CO.dif': [((('E', 1), ('CO', 2)), 1), ((('E', 1), ('CO', 2), ('E', 3)), 2)],
            'k.COO.lh': [],
            'k.COO.el': []
        }

        # Get the dictionary of states.
        decay_dict_resultant1 = system._get_multiplicity(decay_dict0)

    def test_get_numbering(self):
        """ Tests that the _get_numbering function is working properly.
        """

        # Define the states and the number of sites.
        states = {'E', 'CO', 'O'}
        number_of_sites = 6

        # Create a CO oxidation system with six sites.
        system = Co(number_of_sites, states)

        # ----------------------------------------------------------------------
        # Get requested states for a 2 state site.
        # ----------------------------------------------------------------------

        # Get the state to be labeled.
        state = ('CO', 'CO')

        # Get the outcomes.
        outcomes_0 = system._get_numbering(state)

        # These are the states that must come out.
        outcomes_1 = [(('CO', 1), ('CO', 2),), (('CO', 2), ('CO', 3),),
                      (('CO', 3), ('CO', 4),), (('CO', 4), ('CO', 5),),
                      (('CO', 5), ('CO', 6),)
                      ]

        # Check that all the states are unique.
        self.assertEqual(len(outcomes_0), len(set(outcomes_0)))
        self.assertEqual(len(outcomes_1), len(set(outcomes_1)))

        # Check that all the outcomes are in the list.
        for outcome_0 in outcomes_0:
            self.assertTrue(outcome_0 in outcomes_1)

        # ----------------------------------------------------------------------
        # Get requested states for a 3 state site.
        # ----------------------------------------------------------------------

        # Get the state to be labeled.
        state = ('CO', 'CO', 'E')

        # Get the outcomes.
        outcomes_0 = system._get_numbering(state)

        # These are the states that must come out.
        outcomes_1 = [(('CO', 1), ('CO', 2), ('E', 3),),
                      (('CO', 2), ('CO', 3), ('E', 4),),
                      (('CO', 3), ('CO', 4), ('E', 5),),
                      (('CO', 4), ('CO', 5), ('E', 6),)
                      ]

        # Check that all the states are unique.
        self.assertEqual(len(outcomes_0), len(set(outcomes_0)))
        self.assertEqual(len(outcomes_1), len(set(outcomes_1)))

        # Check that all the outcomes are in the list.
        for outcome_0 in outcomes_0:
            self.assertTrue(outcome_0 in outcomes_1)

        # ----------------------------------------------------------------------
        # Get requested states for a 7 state site; not possible.
        # ----------------------------------------------------------------------

        # Get the state to be labeled with an excess of entries.
        state = ('CO', 'CO', 'E', 'E', 'E', 'E', 'E',)

        # The request must fail.
        try:
            # Get the outcomes.
            system._get_numbering(state)
            self.assertTrue(False)
        except ValueError:
            self.assertTrue(True)

        # Get the state to be labeled with zero entries.
        state = []

        # The request must fail.
        try:
            # Get the outcomes.
            system._get_numbering(state)
            self.assertTrue(False)
        except ValueError:
            self.assertTrue(True)

    def test_get_states(self):
        """ Tests that the _get_states function is working properly.
        """

        # Define the states and the number of sites.
        states = {'E', 'CO', 'O'}
        number_of_sites = 6

        # Create a CO oxidation system with six sites.
        system = Co(number_of_sites, states)

        # ----------------------------------------------------------------------
        # Get requested states for a 2 state site.
        # ----------------------------------------------------------------------

        # Define the order.
        order = 2

        # Get the states.
        outcomes_0 = system._get_states(order)

        # These are the states that must come out.
        outcomes_1 = [('CO', 'CO',), ('CO', 'O',), ('CO', 'E',),
                      ('O', 'CO',), ('O', 'O',), ('O', 'E',),
                      ('E', 'CO',), ('E', 'O',), ('E', 'E',)
                      ]

        # Check that all the states are unique.
        self.assertEqual(len(outcomes_0), len(set(outcomes_0)))
        self.assertEqual(len(outcomes_1), len(set(outcomes_1)))

        # Check that all the outcomes are in the list.
        for outcome_0 in outcomes_0:
            self.assertTrue(outcome_0 in outcomes_1)


if __name__ == '__main__':
    # Run the tests.
    unittest.main()

    # Reset the path.
    os.chdir(path1)
    print(os.getcwd())
