# Imports
import copy as cp
import itertools
import os


class EquationGenerator:

    # --------------------------------------------------------------------------
    # Auxiliary methods.
    # --------------------------------------------------------------------------

    @staticmethod
    def get_all_equations(destination_path, file_name="equations"):
        """ Generates all the equations and saves them in a file in LaTeX form
            to copy and paste in the document.

            :param destination_path: The full destination path of where the file
            is to be saved.

            :param file_name: The name of the file where the equations are to be
             saved. The name must be provided without an extension.
        """

        # ----------------------------------------------------------------------
        # File properties setup.
        # ----------------------------------------------------------------------

        # Make sure the path separator works.
        full_path = destination_path.strip()
        full_path += "" if full_path[-1] == os.sep else os.sep

        # Strip the name from spaces and/or extensions.
        name = file_name.strip()
        name = name.split(".")[0]
        name = name.strip()

        # Get the full path.
        full_file_path = full_path + name

        # ----------------------------------------------------------------------
        # Generate and save the quantities.
        # ----------------------------------------------------------------------

        # All the possible states for a site.
        site_states = ["CO", "O", "E"]

        # All the possible states in the system.
        all_states = [[i, j, k] for i, j, k in itertools.product(site_states, site_states, site_states)]

        # Do it for all the states.
        for i, state in enumerate(all_states):
            # Number the file correctly.
            file_path = full_file_path + f"({i+1}).txt"

            # Get all the incoming terms and outgoing terms.
            table_input = EquationGenerator.generate_incoming_rates(state)
            table_outpu = EquationGenerator.generate_outgoing_rates(state)

            # Remember to copy the remaining output to the table.
            table_input[i][0], table_input[i][1] = table_outpu[0][0], table_outpu[0][1]

            # Generate the HTML table.
            table_string = EquationGenerator.generate_table(state, table_input, f"{i + 1}. ", 2)

            # Generate the LaTeX format equation.
            equation_string = EquationGenerator.generate_equation(state, table_input)

            # Merge the table and equation.
            state_string = table_string + "\n\n***\n\n" + equation_string

            # Save the variable to the file.
            with open(file_path, "w") as fl:
                fl.write(state_string + "\n")

    @staticmethod
    def get_latex_form(state):
        """ Given a state, it returns its latex type form.

            :param state: The given valid state.

            :return latex_string: The LaTeX string representation of the state.
        """

        # Add the left bracket.
        latex_string = f"\left<"

        # Add the site states with the appropriate numbering.
        for i, site in enumerate(state):
            latex_string += str(site) + "_{" + str(i + 1) + "}"

        # Add the right bracket.
        latex_string += f"\\right>"

        return latex_string

    @staticmethod
    def get_gain_rate(initial_state, final_state):
        """ Get the total rate at which the initial state generates the final
            state.

            :param initial_state: The initial state from which the final state
            must be generated.

            :param final_state: The final state to which the initial state is
            intended to decay.

            :return: A string representing the total rate at which the initial
            state decays into the final state.
        """

        # Auxiliary variables.
        rates_list = []

        # String where the total rate is stored.
        rates_string = ""

        # States that can be, potentially, generated from the initial state.
        rates, possible_states = EquationGenerator.get_rate_function_tuple(initial_state)

        # Search for the process that generate the final state.
        for i, states in enumerate(possible_states):
            # If final state is NOT in the list, continue.
            if final_state not in states:
                continue

            # Add the rate to the possible outcomes.
            rates_list.append(rates[i])

        # Check that there is at least one rate to continue.
        if len(rates_list) == 0:
            return "0"

        # Generate the string.
        rates_string = rates_list[0] if len(rates_list) == 1 else "\\left(" + " + ".join(rates_list) + "\\right)"

        return rates_string

    @staticmethod
    def get_rate_function_tuple(state):
        """ Returns a tuple of tuples that has the symbolic representation of
            the reaction constant and the list of the resultant configurations
            of applying the given processes to the given state.

            :param state: The state to which the functions are applied.

            :return information_tuple: A tuple that contains the symbolic
            representation of the reaction constant and the lists of the
            resultant configurations of applying the given processes to the
            given state.
        """

        # Define the tuple with the rate symbols of the different processes.
        rate_labels = (
            "k_{O_{ads}}",
            "k_{O_{des}}",
            "k_{O_{dif}}",
            "k_{CO_{ads}}",
            "k_{CO_{des}}",
            "k_{CO_{dif}}",
            "k_{COO_{l-h}}",
            "k_{COO_{e-r}}"
        )

        # Get the lists of the final states generated from each process applied
        # to the given state.
        final_states = (
            EquationGenerator.oxygen_adsorb(state),
            EquationGenerator.oxygen_desorb(state),
            EquationGenerator.oxygen_diffusion(state),
            EquationGenerator.carbon_monoxide_adsorb(state),
            EquationGenerator.carbon_monoxide_desorb(state),
            EquationGenerator.carbon_monoxide_diffusion(state),
            EquationGenerator.lh_reaction(state),
            EquationGenerator.er_reaction(state)
        )

        return rate_labels, final_states

    @staticmethod
    def get_unique(states):
        """ Makes the list of lists such that they have unique elements.

            :param states: The list of states to be made unique.

            :return unique_states: A list with unique elements.
        """

        # Make elements hashable and set them.
        unique_states = set(map(tuple, states))

        # convert all lists into lists again.
        unique_states = list(map(list, unique_states))

        return unique_states

    # --------------------------------------------------------------------------
    # Generate rates.
    # --------------------------------------------------------------------------

    @staticmethod
    def generate_equation(state, table):
        """ Generates a rate equation with the given transition rates in the
            table.

            :param state: The state for which the differential equation is to be
            written.

            :param table: The table that contains the information of the  rates
            and states of the system.

            :return eq_string: A string representing an equation with the
            contribution of each state to the overall gain/decay rate of a
            state.
        """

        # Generate a list with the non-zero terms.
        non_zero_rates = [ "".join(x) for x in table if not x[0] == "0"]

        # Join the strings.
        eq_string = "$$\n\t\\begin{align}\n"
        eq_string += "\t\t\\frac{d" + EquationGenerator.get_latex_form(state) + "}{dt}=&"
        for i, strng in enumerate(non_zero_rates):
            eq_string += " + " + strng if i > 0 and not strng[0] == "-" else strng
        eq_string += "\n\t\\end{align}\n$$"

        return eq_string

    @staticmethod
    def generate_incoming_rates(final_state):
        """ Generates the string that represents the total rate at which the
            given state is generated from other states.

            :param final_state: The state intended to be reached by means of the
            possible process on other states.

            :return rates_string: The string that represents TOTAL sum of the
            states that can generate the given state.
        """

        # Possible states.
        states = ["CO", "O", "E"]

        # Generate a list of all possible states.
        all_states = [[i, j, k] for i, j, k in itertools.product(states, states, states)]

        # List that contains the rate contribution from each state to producing
        # the final state.
        rate_list = [["0", ""] for _ in range(len(all_states))]

        # Scan all states.
        for i, initial_state in enumerate(all_states):

            # Fill in the LaTeX form of the name from the list.
            rate_list[i][1] = EquationGenerator.get_latex_form(initial_state)

            # Do it for states different from the final_state.
            if final_state == initial_state:
                continue

            # Get the total rate at which the initial state generates the final state.
            rate_list[i][0] = EquationGenerator.get_gain_rate(initial_state, final_state)

        return rate_list

    @staticmethod
    def generate_outgoing_rates(initial_state):
        """ Generates the string that represents the total rate at which the
            given state decays.

            :param initial_state: The state from which the decay will be
            generated.

            :return rates_string: The string that represents TOTAL sum of the
            rates at which the given state decays.
        """

        # Auxiliary variables.
        tmp_list = []
        rate_list = [["0", ""]]

        # Get the string representation of the state.
        latex_string = EquationGenerator.get_latex_form(initial_state)
        rate_list[0][1] = latex_string

        # Get the functions to apply.
        labels, final_states = EquationGenerator.get_rate_function_tuple(initial_state)

        # Get the constants into a list.
        for i, states0 in enumerate(final_states):
            # Get the number of states.
            state_number = len(states0)

            # No need to include the lists of zero length.
            if state_number == 0:
                continue

            # Get the proper string.
            tmp = "" if state_number == 1 else f"{state_number}"
            tmp += labels[i]

            # Append the string to the rates string.
            tmp_list.append(tmp)

        # If there are no rates to append return an empty string.
        if len(tmp_list) == 0:
            return rate_list

        # Generate the parenthesized expression.
        rates_string = f"-{tmp_list[0]}" if len(tmp_list) == 1 else "-\\left(" + " + ".join(tmp_list) + "\\right)"

        # Save the rate and return the list.
        rate_list[0][0] = rates_string

        return rate_list

    @staticmethod
    def generate_table(state, table, index="", columns=6):
        """ Generates an HTML table with the given transition rates in the
            table.

            :param state: The state for which the table is to be generated.

            :param table: The table that contains the information of the
             rates and states of the system.

            :param index: The index with which the entry must be labeled.

            :param columns: The number of columns the table should have.

            :return tb_string: A string representing an HTML table with the
            contribution of each state to the overall gain/decay rate of a
            state.
        """

        # Auxiliary variables.
        tb_header = str(index) + "$\\frac{d" + EquationGenerator.get_latex_form(state) + "}{dt}$:\n\n"
        tb_header += '<table style="background-color: white; width: 100%; text-align: left">\n'
        tb_footer = '</table>'

        td_opening = '\t\t<td style="background-color: white; width: 100%; text-align: left">\n'
        td_closing = '\t\t</td>\n'

        tr_opening = '\t<tr>\n'
        tr_closing = '\t</tr>\n'

        # Start with the header and an opening table row.
        tb_string = tb_header + tr_opening

        for i, term in enumerate(table):
            if i % columns == 0 and i > 0:
                tb_string += tr_closing + tr_opening

            tb_string += td_opening
            tb_string += "\t\t\t" + str(i + 1) + ". $" + term[0] + "$\n"
            tb_string += td_closing

            # Complete the table.
            if i == len(table) - 1:
                for j in range((i+1) % columns):
                    tb_string += td_opening + td_closing

        # Remember to always close the table.
        tb_string += tr_closing + tb_footer

        return tb_string

    # --------------------------------------------------------------------------
    # Oxygen exclusive methods.
    # --------------------------------------------------------------------------

    @staticmethod
    def oxygen_adsorb(initial_state):
        """ Given a state, returns the list of states that are generated by the
            adsorption of oxygen.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of oxygen adsorption.
        """

        # Determine if the state is valid for analysis.
        EquationGenerator.validate_state(initial_state)

        # Define the list of possible final states.
        final_states = []

        # Check both sites
        for i in range(0, 2):
            # Auxiliary variables.
            j = i + 1

            # Always make a copy of the list.
            tmp = cp.deepcopy(initial_state)

            # Check which sites can adsorb oxygen and do the process.
            if tmp[i] == "E" and tmp[j] == "E":
                tmp[i], tmp[j] = "O", "O"
                final_states.extend([tmp])

        # Get the unique elements.
        final_states = EquationGenerator.get_unique(final_states)

        # Determine if the final states are valid.
        for state in final_states:
            EquationGenerator.validate_state(state)

        return final_states

    @staticmethod
    def oxygen_desorb(initial_state):
        """ Given a state, returns the list of states that are generated by the
            desorption of oxygen.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of oxygen desorption.
        """

        # Determine if the state is valid for analysis.
        EquationGenerator.validate_state(initial_state)

        # Define the list of possible final states.
        final_states = []

        # Check both sites
        for i in range(0, 2):
            # Auxiliary variables.
            j = i + 1
            # Always make a copy of the list.
            tmp = cp.deepcopy(initial_state)

            # Check which sites can desorb oxygen and do the process.
            if tmp[i] == "O" and tmp[j] == "O":
                tmp[i], tmp[j] = "E", "E"
                final_states.extend([tmp])

        # Get the unique elements.
        final_states = EquationGenerator.get_unique(final_states)

        # Determine if the final states are valid.
        for state in final_states:
            EquationGenerator.validate_state(state)

        return final_states

    @staticmethod
    def oxygen_diffusion(initial_state):
        """ Given a state, returns the list of states that are generated by the
            diffusion of oxygen to adjacent sites.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of oxygen diffusion to adjacent sites.
        """

        # Determine if the state is valid for analysis.
        EquationGenerator.validate_state(initial_state)

        # Define the list of possible final states.
        final_states = []

        # Check both sites
        for i in range(0, 2):
            # Auxiliary variables.
            j = i + 1
            # Always make a copy of the initial state.
            tmp = cp.deepcopy(initial_state)

            # Try to make the diffusion happen.
            if (tmp[i] == "O" and tmp[j] == "E") or (tmp[i] == "E" and tmp[j] == "O"):
                # Exchange the particles.
                tmp[i], tmp[j] = tmp[j], tmp[i]

                # Append to the list.
                final_states.extend([tmp])

        # Get the unique elements.
        final_states = EquationGenerator.get_unique(final_states)

        # Determine if the final states are valid.
        for state in final_states:
            EquationGenerator.validate_state(state)

        return final_states

    # --------------------------------------------------------------------------
    # Carbon monoxide exclusive methods.
    # --------------------------------------------------------------------------

    @staticmethod
    def carbon_monoxide_adsorb(initial_state):
        """
            Given a state, returns the list of states that are generated by the
            adsorption of carbon monoxide.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of carbon monoxide adsorption.
        """

        # Determine if the state is valid for analysis.
        EquationGenerator.validate_state(initial_state)

        # Define the list of possible final states.
        final_states = []

        # Check ALL sites
        for i in range(0, 3):
            # Always make a copy of the list.
            tmp = cp.deepcopy(initial_state)

            # Check which sites can adsorb carbon monoxide and do the process.
            if tmp[i] == "E":
                tmp[i] = "CO"
                final_states.extend([tmp])

        # Get the unique elements.
        final_states = EquationGenerator.get_unique(final_states)

        # Determine if the final states are valid.
        for state in final_states:
            EquationGenerator.validate_state(state)

        return final_states

    @staticmethod
    def carbon_monoxide_desorb(initial_state):
        """
            Given a state, returns the list of states that are generated by the
            desorption of carbon monoxide.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of carbon monoxide desorption.
        """

        # Determine if the state is valid for analysis.
        EquationGenerator.validate_state(initial_state)

        # Define the list of possible final states.
        final_states = []

        # Check both sites
        for i in range(0, 3):
            # Always make a copy of the list.
            tmp = cp.deepcopy(initial_state)

            # Check which sites can desorb carbon monoxide and do the process.
            if tmp[i] == "CO":
                tmp[i] = "E"
                final_states.extend([tmp])

        # Get the unique elements.
        final_states = EquationGenerator.get_unique(final_states)

        # Determine if the final states are valid.
        for state in final_states:
            EquationGenerator.validate_state(state)

        return final_states

    @staticmethod
    def carbon_monoxide_diffusion(initial_state):
        """
            Given a state, returns the list of states that are generated by the
            diffusion of carbon monoxide to adjacent sites.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the process of carbon monoxide diffusion to adjacent sites.
        """

        # Determine if the state is valid for analysis.
        EquationGenerator.validate_state(initial_state)

        # Define the list of possible final states.
        final_states = []

        # Check both sites
        for i in range(0, 2):
            # Auxiliary variables.
            j = i + 1

            # Always make a copy of the initial state.
            tmp = cp.deepcopy(initial_state)

            # Try to make the diffusion happen.
            if (tmp[i] == "CO" and tmp[j] == "E") or (tmp[i] == "E" and tmp[j] == "CO"):
                # Exchange the particles.
                tmp[i], tmp[j] = tmp[j], tmp[i]

                # Append to the list.
                final_states.extend([tmp])

        # Get the unique elements.
        final_states = EquationGenerator.get_unique(final_states)

        # Determine if the final states are valid.
        for state in final_states:
            EquationGenerator.validate_state(state)

        return final_states

    # --------------------------------------------------------------------------
    # Carbon monoxide - oxygen reaction exclusive methods.
    # --------------------------------------------------------------------------

    @staticmethod
    def lh_reaction(initial_state):
        """
            Given a state, returns the list of states that are generated by the
            Langmuir-Hinshelwood reaction of carbon monoxide and oxygen; i.e.,
            oxygen and carbon monoxide on adjacent sites react and immediately
            desorb.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to the Langmuir-Hinshelwood reaction of carbon monoxide and oxygen.
        """

        # Determine if the state is valid for analysis.
        EquationGenerator.validate_state(initial_state)

        # Define the list of possible final states.
        final_states = []

        # Check ALL sites
        for i in range(0, 2):
            # Auxiliary variables.
            j = i + 1
            # Always make a copy of the list.
            tmp = cp.deepcopy(initial_state)

            # Try to make the reaction happen.
            if (tmp[i] == "CO" and tmp[j] == "O") or (tmp[i] == "O" and tmp[j] == "CO"):
                tmp[i], tmp[j] = "E", "E"
                final_states.extend([tmp])

        # Get the unique elements.
        final_states = EquationGenerator.get_unique(final_states)

        # Determine if the final states are valid.
        for state in final_states:
            EquationGenerator.validate_state(state)

        return final_states

    @staticmethod
    def er_reaction(initial_state):
        """
            Given a state, returns the list of states that are generated by the
            gas-phase reaction oxygen on the surface with gaseous carbon
            monoxide in the environment; i.e., a surface oxygen and gaseous
            carbon monoxide Elay-Rideal reaction.

            :param initial_state: The state to be analyzed.

            :return: final_states: A list of final states that are generated due
            to a surface oxygen and gaseous carbon monoxide Elay-Rideal
            reaction.
        """

        # Determine if the state is valid for analysis.
        EquationGenerator.validate_state(initial_state)

        # Define the list of possible final states.
        final_states = []

        # Check ALL sites
        for i in range(0, 3):
            # Always make a copy of the list.
            tmp = cp.deepcopy(initial_state)

            # Check which sites can have an Elay-Rideal reaction.
            if tmp[i] == "O":
                tmp[i] = "E"
                final_states.extend([tmp])

        # Get the unique elements.
        final_states = EquationGenerator.get_unique(final_states)

        # Determine if the final states are valid.
        for state in final_states:
            EquationGenerator.validate_state(state)

        return final_states

    # --------------------------------------------------------------------------
    # Validation methods.
    # --------------------------------------------------------------------------

    @staticmethod
    def validate_state(state):
        """ Validates that the state is a list of three items and those items
            are in the list ["CO", "O", "E"].
        """

        # List of possible states.
        state_list = ["CO", "O", "E"]

        # The state to be validated.
        if not isinstance(state, (list,)):
            raise ValueError(f"The type of the state is not valid. Current  type: {type(state)}.")

        elif not len(state) == 3:
            raise ValueError(f"The current length of the state lis is not valid, it must be 3."
                             f" Current legth: {len(state)}.")

        elif not all(map(lambda x: x in state_list, state)):
            raise ValueError(f"The states are not valid, they must be in the list {state_list}."
                             f" Current state: {state}.")


if __name__ == "__main__":

    CD = os.path.dirname(__file__)

    EquationGenerator.get_all_equations(CD)
