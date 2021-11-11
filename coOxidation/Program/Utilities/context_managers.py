""" File that contains the different context managers."""

# ------------------------------------------------------------------------------
# Imports.
# ------------------------------------------------------------------------------

# Imports: General.
import copy
import os

# ------------------------------------------------------------------------------
# Classes.
# ------------------------------------------------------------------------------


class WorkingDirectorySet:
    """ Context manager that sets the current working directory (CWD).

        - self.oldpath: Variable that stores the old path.

        - self.newpath: Variable that stores the new path.
    """

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Constructors and Dunder Methods.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # --------------------------------------------------------------------------
    # Constructor(s).
    # --------------------------------------------------------------------------

    def __init__(self, newpath: str):
        """ Initializes the context manager.

            :param newpath: The path of the directory to be set.
        """

        self.oldpath = copy.deepcopy(os.getcwd())
        self.newpath = newpath

    # --------------------------------------------------------------------------
    # Dunder Methods.
    # --------------------------------------------------------------------------

    def __enter__(self) -> str:
        """ Sets the working directory and returns the directory path.

            :return: The directory path.
        """

        self.newpath = self.newpath.strip()
        self.newpath = self.newpath if self.newpath[-1] == os.sep else f"{self.newpath}{os.sep}"
        os.chdir(self.newpath)
        return self.newpath

    def __exit__(self, exc_type, exc_val, exc_tb):
        """ Sets the working directory to the initial path.
        """
        os.chdir(self.oldpath)
