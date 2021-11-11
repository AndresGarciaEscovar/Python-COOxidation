import numpy

from coOxidation.Program.KMC.COOxidation_parameters import COOxidationKMCParameters
from coOxidation.Program.KMC.COOxidation_KMC import COOxidationKMC

if __name__ == "__main__":

    parameters = COOxidationKMCParameters()
    parameters.length = 3
    parameters.maximum_counter = 9

    simulation = COOxidationKMC(parameters)
    simulation.reset_simulation(6.0)

