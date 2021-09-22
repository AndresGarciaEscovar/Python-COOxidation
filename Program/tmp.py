
# Imports: Class to be tested.
from Program.co_oxidation_model import COOxidationEquationGenerator as Co

if __name__ == "__main__":
    # Define the states and the number of sites.
    states = {'E', 'CO', 'O'}
    number_of_sites = 3

    # Create a CO oxidation system with three sites.
    system = Co(number_of_sites, states)

    # system.get_nth_order_equations(system.number_of_sites)
    operations = system._get_associated_operations()
    system.get_nth_order_equations(operations, order=2)
