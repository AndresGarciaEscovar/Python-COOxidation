""" File that contains the class to define the simulation parameters."""

# ------------------------------------------------------------------------------
# Imports.
# ------------------------------------------------------------------------------

# Imports: General.
import csv

from matplotlib import pyplot

# ------------------------------------------------------------------------------
# Classes.
# ------------------------------------------------------------------------------


class COOxidationAnalysis:
    """ Class that contains the analysis tools for the simulation.
    """

    @staticmethod
    def plot_results(file_path: str) -> None:
        """ Plots the results of the simulation.

            :param file_path: The path of the file where the results are saved.
        """

        with open(file_path, "r", newline="\n") as fl:
            reader = csv.reader(fl, delimiter=",")
            data = [row for row in reader]

        info = data[:2]
        data = data[2:-1]

        labels = tuple(label for label in info[1][1:])
        fig, axes = pyplot.subplots(ncols=len(data))
        for i, axis in enumerate(data):
            data_ = data[i][1:]
            axes[i].title.set_text(f"Site {i + 1}")
            axes[i].pie(data_, labels=labels, autopct='%1.5f%%', normalize=True)

        n = 6
        title = [info[0][k: k + n] for k in range(0, len(info[0]), n)]
        title = list(map(lambda x: ", ".join(x), title))
        title = "\n".join(title)

        pyplot.suptitle(title, fontsize=10)
        pyplot.tight_layout()
        pyplot.show()







