import itertools


state = ((('CO', 1),('CO', 2), ('CO', 3)),
         (('CO', 1),('CO', 2), ('O', 3)),
         (('CO', 1),('CO', 2), ('E', 3)),
         (('CO', 1),('O', 2), ('E', 3)),
         )

print((("CO", 1),("CO", 2),("CO", 3)) in state)
