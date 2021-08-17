# Imports: General.
import numpy as np

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

    # Oxygen adsorption.
    k[0] = 1

    # Oxygen desorption.
    k[1] = 1

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


def choose_random_number(number1, number2):
    """ Gets a random number in the range (number1, number2)
    
        :param number1: The first number in the sequence.
        
        :param number2: The second number in the sequence.
        
        :return random_number: A random number in the range (number1, number2)
    """
    import random

    # Choose a random number that does NOT include either of the ends.
    random_number = random.uniform(number1, number2)
    while random_number == number1 or random_number == 2:
        random_number = random.uniform(number1, number2)

    return random_number


def choose_random_move():
    """ Chooses a random move to be performed.
    """

    # Choose a random site.
    site = choose_random_site()

    # Choose a random move, i.e., a random rate from the list.
    rates = get_rates()
    random_number = choose_random_number(0, rates[-1])

    # Check that everything is consistent.
    if random_number > rates[-1]:
        raise ValueError(f"The number chosen is greater than that of the total rate; this should not happen.\n"
                         f"Maximum rate: {rates[-1]}, Chosen rate: {random_number}.")

    # Choose the move.
    for i, rate in rates:
        if random_number < rate:
            rate_id = i
            break

    # --------------------------------------------------------------------------
    # Perform an adsorption/desorption move.
    # --------------------------------------------------------------------------

    # Choose the move.
    if rate_id == 0 or  rate_id == 3: # Adsorption.
        adsorb_desorb_particle(site, rate_id, "adsorb")

    elif rate_id == 1 or rate_id == 4: # Desorption.
        adsorb_desorb_particle(site, rate_id, "desorb")

    elif rate_id == 2 or rate_id == 5: # Diffusion.
        move_particle(site, rate_id)

    elif rate_id == 6: # CO-O reaction on surface.
        reaction_particle(site, rate_id)

    elif rate_id == 7: # CO-O reaction in gas.
        reaction_particle(site, rate_id)

    elif rate_id == 8: # CO-O reaction anywhere.
        reaction_particle(site, rate_id)


if __name__ == "__main__":

    # Create the array where the particles will live.
    particle_array = ["" for _ in range(3)]

    final_state_count = [0 for _ in range(3)]

    # Particle types.
    particle_types = ['E', "O", "CO"]

    choose_random_move()


