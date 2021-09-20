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

a = ['O', 'E', 'CO']

b = ['O', 'CO', 'E']

print(np.array_equal(a,b))


