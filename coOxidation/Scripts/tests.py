from itertools import product


if __name__ == "__main__":
    # Set a fictitious state.
    states = ['CO', 'O', 'E']

    states = product(*([states] * 3))

    for i, state in enumerate(states):
        print(i + 1, ". ", state)

    # TODO: CONTINUE IN THE equation_generator.py FILE.
