import os.path

from coOxidation.Program.KMC.COOxidationAnalysis import COOxidationAnalysis
from coOxidation.Program.KMC.COOxidation_parameters import COOxidationKMCParameters
from coOxidation.Program.KMC.COOxidation_KMC import COOxidationKMC

from coOxidation.Program.Utilities.context_managers import WorkingDirectorySet

if __name__ == "__main__":
    with WorkingDirectorySet(os.path.dirname(__file__)) as fl:
        parameters = COOxidationKMCParameters()
        parameters.length = 3
        parameters.maximum_counter = 0.5
        parameters.repetitions = 100_001

        simulation = COOxidationKMC(parameters)
        simulation.run_simulation()
        simulation.statistics_save("results.txt")

        COOxidationAnalysis.plot_results("results.txt")
