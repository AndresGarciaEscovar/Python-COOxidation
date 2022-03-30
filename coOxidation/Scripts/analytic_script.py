
# Imports: Class to be tested.
import os.path

from coOxidation.Program.Analytic.equation_generator import EquationGenerator

if __name__ == "__main__":
    # Define the states and the number of sites.
    sites_number = 3

    # Define the equation order (i.e., 0 = Exact, 1 = MF, 2 = PA, ...).
    order = 2

    # Create a CO oxidation system with three sites.
    system = EquationGenerator(sites_number)

    # Get the nth order equations .
    system.order = order
    system.get_nth_order_equations(display=True)

    # Generate the equations for the given format.
    system.save_equations(file_name="tmp", format_type="Mathematica", order=order, save_path=os.path.dirname(__file__))
