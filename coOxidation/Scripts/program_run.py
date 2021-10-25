
# Imports: Class to be tested.
import os.path

from coOxidation.Program.Analytic.equation_generator import EquationGenerator

if __name__ == "__main__":
    # Define the states and the number of sites.
    number_of_sites = 3

    # Define the equation order.
    order = 2

    # Create a CO oxidation system with three sites.
    system = EquationGenerator(number_of_sites)

    # Get the nth order equations.
    system.get_nth_order_equations(order=order, print_equations=False)

    # Generate the equations.
    system.save_equations(file_name="tmp", format_type="Mathematica", order=order, save_path=os.path.dirname(__file__))
