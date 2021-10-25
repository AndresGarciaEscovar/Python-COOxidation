from coOxidation.Program.Analytic.Formatters.Formatters.formatter_mathematica import MathematicaFormatter

if __name__ == "__main__":
    # Set a fictitious state.
    state = (("CO", 1), ("CO", 2),)

    print(MathematicaFormatter.get_state(state, 3))

