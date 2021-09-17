import itertools

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
for comb in combs:
    # Only look a the tuples that have elements of the same length.
    if len(set(tuple(map(len, comb)))) == 1:
        candidates.append(comb)

# No more need for this variable.
del combs

for candidate in candidates:
    for elem in candidate:
        print(elem)
    print("")