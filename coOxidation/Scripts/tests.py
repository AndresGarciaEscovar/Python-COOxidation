from coOxidation.Program.Analytic.equation_generator import EquationGenerator

if __name__ == "__main__":
    # Set a fictitious state.
    states = [
        (("CO", 1), ("CO", 2), ("CO", 3),),
        (("CO", 1), ("CO", 2), ("E", 3),),
        (("CO", 1), ("CO", 2), ("O", 3), )
    ]

    print(EquationGenerator(3).get_nth_order_equations(True))



