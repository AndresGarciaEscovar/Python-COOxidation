
# Imports: Class to be tested.
import os.path

from Program.co_oxidation_model import COOxidationEquationGenerator as Co

if __name__ == "__main__":
    # Define the states and the number of sites.
    states = {'E', 'CO', 'O'}
    number_of_sites = 3

    # Create a CO oxidation system with three sites.
    system = Co(number_of_sites, states)

    # Get the nth order equations.
    system.get_nth_order_equations(order=3, print_equations=False)

    # Generate the equations.
    system.get_equations_in_format(file_name="tmp", format_type="latex", order=0, save_path=os.path.dirname(__file__), together=True)
