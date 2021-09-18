import itertools
import numpy as np

# Define the possible states.
possible_states = sorted(['E', 'O', 'CO'])

state_list = [(('E', 2), ('E', 3)),
              (('E', 2), ('E', 3)),
              (('CO', 1), ('E', 2), ('E', 3)),
              (('CO', 1), ('E', 2), ('E', 3)),
              (('O', 1), ('E', 2), ('E', 3)),
              (('O', 1), ('E', 2), ('E', 3)),
              (('E', 1), ('E', 2), ('E', 3)),
              (('E', 1), ('E', 2), ('E', 3)),
              (('E', 1), ('E', 2), ('E', 3))]

# Auxiliary variables.
candidates = []

# Gather the tuples in groups of three.
combs = list(itertools.combinations(state_list, 3))

# Search the possible combinations.
for comb in combs:
    # Get the lengths of the all the states in the combination.
    length_list = (len(x) for x in comb)

    # Only look a the tuples that have elements of the same length.
    if len(set(length_list)) == 1:
        candidates.append(comb)

# No more need for this variable anymore.
del combs

# Candidate entry.
candidate_entry = -1
reduced_site = []

# Get the indexes and the states for each tuple.
for i, candidate_states in enumerate(candidates):
    # Auxiliary variables.
    sites_indexes = []
    sites_states = []

    # Get the indexes of each tuple.
    for state in candidate_states:
        # Separate the indexes and the states.
        indexes = tuple(x[1] for x in state)
        states = tuple(x[0] for x in state)

        # Append them to the lists.
        sites_indexes.append(indexes)
        sites_states.append(states)

    # Only states that have the same indexes can be examined.
    if len(set(sites_indexes)) == 1:
        # Convert into a numpy array and analyze the states.
        sites_states = np.array(sites_states)

        # Analyze the ends to see if the state can be reduced.
        if sorted(list(sites_states[:, 0])) == possible_states:
            candidate_entry = i
            reduced_site = tuple(x for j, x in enumerate(candidates[candidate_entry][0]) if j > 0)
            break

# Get the indexes of the states that are to be removed.
index_list = []
for state in candidates[candidate_entry]:
    index_list.append(state_list.index(state))

print(index_list)


for index in reversed(index_list):
    state_list.pop(index)


for state in state_list:
    print(state)



