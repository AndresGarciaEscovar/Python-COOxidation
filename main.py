# Imports: General.
import random

import numpy as np

# ------------------------------------------------------------------------------
# Global Variables.
# ------------------------------------------------------------------------------

# The  configuration array of the system
configuration_array = ['E', 'E', 'E']

# A dictionary with the particle types and their final count.
final_state_count = [{"E": 0, "O": 0,"CO": 0} for _ in range(3)]

# Elapsed time.
elapsed_time = np.double(0.0)

# Maximum time, i.e., time to run the simulation.
maximum_time = np.double(10 ** 3)

# Number of simulations over which to average.
step_number = 10

# ------------------------------------------------------------------------------
# Get Functions.
# ------------------------------------------------------------------------------


def get_rates():
    """ Gets the rates of the system. More rates can be added if needed.

        :return k: Returns the cumulative rates of the system.
    """

    # There are 9 rates in total.
    k = [0 for _ in range(9)]

    # --------------------------------------------------------------------------
    # Oxygen related rates.
    # --------------------------------------------------------------------------

    # Oxygen adsorption, always try to adsorb on neighboring sites (left or
    # right) with a rate of 1.
    k[0] = 2 * 1

    # Oxygen desorption, always try to desorb from neighboring sites (left or
    # right) with a rate of 1.
    k[1] = 2 * 1

    # Oxygen diffusion, two possible moves (left or right) with a rate of 1.
    k[2] = 2 * 1

    # --------------------------------------------------------------------------
    # Carbon monoxide related rates.
    # --------------------------------------------------------------------------

    # Carbon monoxide adsorption.
    k[3] = 1

    # Carbon monoxide desorption.
    k[4] = 1

    # Carbon monoxide diffusion, two possible moves (left or right) with a rate of 1.
    k[5] = 2 * 1

    # --------------------------------------------------------------------------
    # Carbon monoxide - oxygen interaction related rates.
    # --------------------------------------------------------------------------

    # Carbon monoxide - oxygen neighbor reaction on surface, two possible
    # reaction sites (left or right) with a rate of 1.
    k[6] = 2 * 1

    # Carbon monoxide - oxygen reaction in gas (i.e., oxygen on surface binds to
    # CO in gas)
    k[7] = 1

    # Carbon monoxide - oxygen reaction anywhere in lattice, two possible
    # reaction sites (i.e., not the chosen site) with a rate of 1.
    k[8] = 2 * 1

    # --------------------------------------------------------------------------
    # Calculate the cumulative rates.
    # --------------------------------------------------------------------------

    k = np.array([sum(k[0:i + 1]) if i > 0 else k[i] for i, _ in enumerate(k)], dtype=np.double)

    return k

# ------------------------------------------------------------------------------
# System move functions.
# ------------------------------------------------------------------------------


def adsorb_desorb_particle(site_0, rate_id_0):
    """ Tries to adsorb a carbon monoxide atom or an oxygen pair, or desorbs a
        an oxygen atom or a carbon monoxide atom.

        :param site_0: The site at which the process will take place.

        :param rate_id_0: The rate that was chosen; determines the process
        adsorption or desorption and the particle type that will undergo the
        process.
    """

    global configuration_array

    # Check for adsortion/desorption of oxygen.
    if rate_id_0 == 0 or rate_id_0 == 1:

        # Choose a neighbor site.
        site_1 = site_0 + 1 if choose_random_number(0, 1) <= 0.5 else site_0 - 1

        # Check that the site is NOT out of bounds.
        if site_1 < 0 or site_1 > 2:
            return

        if rate_id_0 == 0 and configuration_array[site_0] == "E" and configuration_array[site_1] == "E":
            configuration_array[site_0] = "O"
            configuration_array[site_1] = "O"
            return

        if rate_id_0 == 1 and configuration_array[site_0] == "O" and configuration_array[site_1] == "O":
            configuration_array[site_0] = "E"
            configuration_array[site_1] = "E"
            return

    # Check for adsortion/desorption of carbon monoxide.
    elif rate_id_0 == 3 or rate_id_0 == 4:

        if rate_id_0 == 3 and configuration_array[site_0] == "E":
            configuration_array[site_0] = "CO"
            return

        if rate_id_0 == 4 and configuration_array[site_0] == "CO":
            configuration_array[site_0] = "E"
            return


def move_particle(site_0, rate_id_0):
    """ Tries to hop a carbon monoxide atom, or an oxygen atom, to a nearest
        neighbor side with equal probability.

        :param site_0: The site at which the process will take place.

        :param rate_id_0: The rate that was chosen; determines which atom is the
        one that will be hopping.
    """

    global configuration_array

    # Choose a random direction.
    site_1 = site_0 + 1 if choose_random_number(0, 1) <= 0.5 else site_0 - 1

    # Check that the site is NOT out of bounds.
    if site_1 < 0 or site_1 > 2:
        return

    # Conditions under which a swap is possible.
    cond_0 = rate_id_0 == 2 and configuration_array[site_0] == "O" and configuration_array[site_1] == "E"
    cond_1 = rate_id_0 == 5 and configuration_array[site_0] == "CO" and configuration_array[site_1] == "E"

    # Make the particle exchange.
    if cond_0 or cond_1:
        configuration_array[site_0], configuration_array[site_1] = configuration_array[site_1], configuration_array[site_0]


def reaction_particle(site_0, rate_id_0):
    """ Reacts a particle or a pair of particles according to the chosen rate,
        i.e., empties the site(s).

        :param site_0: The site at which the process will take place.

        :param rate_id_0: The rate that was chosen; determines which desorption
        process will be attempted.
    """

    global configuration_array

    # --------------------------------------------------------------------------
    # Auxiliary Functions.
    # --------------------------------------------------------------------------

    def validate_sites(site_1_0, rate_id_1_0):
        """ Validates that the site and reaction index are in the proper range.

            :param site_1_0: The site at which the process will take place.

            :param rate_id_1_0: The rate that was chosen; determines which desorption
            process will be attempted.
        """

        # Validate the new site.
        if site_1_0 not in [0, 1, 2]:
            raise ValueError(f"The requested site for reaction is not valid, it must have a value between 0 and 2."
                             f"site_0: {site_1_0}.")

        # Validate the reaction index.
        if rate_id_1_0 not in [6, 7, 8]:
            raise ValueError(f"The requested reaction index is not valid, it must have a value between 6 and 8."
                             f"rate_id_0: {rate_id_1_0}.")

    # --------------------------------------------------------------------------
    # Implementation.
    # --------------------------------------------------------------------------

    # Validate the indexes.
    validate_sites(site_0, rate_id_0)

    # Carbon monoxide - oxygen neighbor reaction on surface, two possible
    # reaction sites (left or right) with a rate of 1.
    if rate_id_0 == 6:
        # Choose a neighboring site.
        site_1 = site_0 + 1 if choose_random_number(0, 1) <= 0.5 else site_0 - 1

        # Check that the site is NOT out of bounds.
        if site_1 < 0 or site_1 > 2:
            return

        # Conditions for desorption.
        cond_0 = configuration_array[site_0] == "CO" and configuration_array[site_1] == "O"
        cond_1 = configuration_array[site_0] == "O" and configuration_array[site_1] == "CO"

        # Desorb if needed.
        if cond_0 or cond_1:
            configuration_array[site_0] = "E"
            configuration_array[site_1] = "E"

        return

    # Carbon monoxide - oxygen reaction in gas (i.e., oxygen on surface binds to
    # CO in gas)
    if rate_id_0 == 7:
        # Desorb the oxygen.
        if configuration_array[site_0] == "O":
            configuration_array[site_0] = "E"

        return

    # Carbon monoxide - oxygen reaction anywhere in lattice, two possible
    # reaction sites (i.e., not the chosen site) with a rate of 1.
    if rate_id_0 == 8:
        # Select the two possible sites to choose.
        sites_0 = [i for i in range(3)]
        sites_0.remove(site_0)

        # Pick a random site.
        site_1 = random.choice(sites_0)

        # Validate the new site.
        if site_1 < 0 or site_1 > 2 or site_1 == site_0:
            raise ValueError(f"The extra site must be 0, 1 or 2 and cannot be the same as site0."
                             f"site_0: {site_0}, site_1: {site_1}")

        # Conditions for desorption.
        cond_0 = configuration_array[site_0] == "CO" and configuration_array[site_1] == "O"
        cond_1 = configuration_array[site_0] == "O" and configuration_array[site_1] == "CO"

        # Desorb if possible.
        if cond_0 or cond_1:
            configuration_array[site_0] = "E"
            configuration_array[site_1] = "E"

        return

# ------------------------------------------------------------------------------
# Random Functions.
# ------------------------------------------------------------------------------


def choose_random_site():
    """
        Chooses an integer between 0 and 2.

        :return: An integer number between 0 and 2.
    """
    import random

    return random.choice([int(x) for x in range(3)])


def choose_random_number(lower_0, upper_0):
    """ Gets a random number in the range (lower_0, upper_0)
    
        :param lower_0: The lower number in the range.
        
        :param upper_0: The upper number in the range.
        
        :return random_number: A random number in the range (lower_0, upper_0)
    """
    import random

    # Choose a random number that does NOT include either of the ends.
    range_list = [lower_0, upper_0]
    random_number = random.uniform(min(range_list), max(range_list))
    while random_number == lower_0 or random_number == upper_0:
        random_number = random.uniform(lower_0, upper_0)

    return random_number


def choose_random_move():
    """ Chooses a random move to be performed.
    """

    # Choose a random site.
    site_0 = choose_random_site()

    # Choose a random move, i.e., a random rate from the list.
    rates_0 = get_rates()
    random_number_0 = choose_random_number(0, rates_0[-1])

    # Check that everything is consistent.
    if random_number_0 > rates_0[-1]:
        raise ValueError(f"The number chosen is greater than that of the total rate; this should not happen.\n"
                         f"Maximum rate: {rates_0[-1]}, Chosen rate: {random_number_0}.")

    # Choose the move.
    rate_id_0 = np.inf
    for i_0, rate_0 in enumerate(rates_0):
        if random_number_0 < rate_0:
            rate_id_0 = i_0
            break

    # --------------------------------------------------------------------------
    # Perform an adsorption/desorption move.
    # --------------------------------------------------------------------------

    # Choose the move.
    if rate_id_0 == 0 or rate_id_0 == 3:  # Adsorption.
        adsorb_desorb_particle(site_0, rate_id_0)

    elif rate_id_0 == 1 or rate_id_0 == 4:  # Desorption.
        adsorb_desorb_particle(site_0, rate_id_0)

    elif rate_id_0 == 2 or rate_id_0 == 5:  # Diffusion.
        move_particle(site_0, rate_id_0)

    elif rate_id_0 == 6:  # CO-O reaction on surface.
        reaction_particle(site_0, rate_id_0)

    elif rate_id_0 == 7:  # CO-O reaction in gas.
        reaction_particle(site_0, rate_id_0)

    elif rate_id_0 == 8:  # CO-O reaction anywhere.
        reaction_particle(site_0, rate_id_0)

    else:
        raise ValueError(f"The chosen rate must be between 0 and 8. Current rate_id: {rate_id_0}")


def increase_time():
    """ Increases the elapsed time by a random amount.
    """

    global elapsed_time

    elapsed_time += -np.log(choose_random_number(0, 1)) / (3 * get_rates()[-1])

    return elapsed_time


if __name__ == "__main__":

    for i in range(step_number):
        # Print a message to the user.
        print(f"Starting simulation {i + 1} of {step_number}:")

        # Always start with the empty configuration array.
        configuration_array = ['E', 'E', 'E']

        # Reset the elapsed time.
        elapsed_time = np.double(0.0)

        # Run the simulation.
        while elapsed_time < maximum_time:
            choose_random_move()
            increase_time()

        # Take the statistics.
        for j in range(3):
            final_state_count[j][configuration_array[j]] += 1

        # Message to the user.
        print(f"Done with simulation {i + 1}.\n")

    # Take the averages.
    for j in range(3):
        for k in ["E", "O", "CO"]:
            final_state_count[j][k] /= step_number

    # Write the results to a file.
    with open("results.txt", "w") as fl:
        fl.write(f"maximum time: {maximum_time}\n")
        for i in [0, 1 , 2]:
            fl.write(f"site {i}:")
            probabilities = [f"P({x}) = {final_state_count[i][x]}" for x in ["E", "O", "CO"]]
            fl.write(", ".join(probabilities) + "\n")



