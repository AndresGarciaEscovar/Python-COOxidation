
class EquationGenerator:

    # --------------------------------------------------------------------------
    # Oxygen exclusive methods.
    # --------------------------------------------------------------------------

    @staticmethod
    def adsorb_oxygen(self, state):

        validate_state(state)

    # --------------------------------------------------------------------------
    # Carbon monoxide exclusive methods.
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # Carbon monoxide - oxygen reaction exclusive methods.
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # Validation methods.
    # --------------------------------------------------------------------------

    @staticmethod
    def validate_state(self, state):
        """ Validates that the state is a list of three items and those items
            are in the list ["E", "O", "CO"].
        """

        if is_instance(state, (list,)) and len(state) == 3:
            print("here")

    # --------------------------------------------------------------------------
    # Constructor.
    # --------------------------------------------------------------------------

    def __init__(self):
        self.tmp = 0


if __name__ == "__main__":

    EquationGenerator.validate_state(["E", "O", "CO"])
