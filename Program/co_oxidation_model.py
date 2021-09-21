""" Writes the equations for carbon monoxide - oxygen associative reaction on
    ruthenium (111).
"""

# Imports: General.
import copy as cp

# Imports: Inherited and auxiliary user-defined classes.
from .mathematica_generator import EquationGenerator

class COOxidation(EquationGenerator):
    """ Writes the equations in various formats, to different orders for the
        carbon monoxide - oxygen associative reaction on ruthenium (111);
        J. Chem. Phys. 143, 204702 (2015). https://doi.org/10.1063/1.4936354

        Class parameters:

        :param self.k_coo_er: String that represents the rate constants of
        carbon monoxide - oxygen gaseous associative desorption of surface
        oxygen; i.e., Elay-Rideal reaction

        :param self.k_coo_lh: String that represents the rate constants of
        carbon monoxide - oxygen neighboring pair associative desorption;
        i.e., Langmuir-Hinshelwoodd reaction.

        :param self.k_co_ads: String that represents the rate constants of
        adsorption of carbon monoxide.

        :param self.k_co_des: String that represents the rate constants of
        desorption of carbon monoxide.

        :param self.k_co_dif: String that represents the rate constants of
        diffusion of carbon monoxide.

        :param self.k_o_ads: String that represents the rate constants of
        adsorption of oxygen.

        :param self.k_o_des: String that represents the rate constants of
        desorption of oxygen.

        :param self.k_o_dif: String that represents the rate constants of
        diffusion of oxygen.

        :param self.o_coo_er: An integer that represents the minimum length of
        the state for carbon monoxide - oxygen gaseous associative desorption of
        surface oxygen; i.e., Elay-Rideal reaction.

        :param self.o_coo_lh: An integer that represents the minimum length of
        the state for carbon monoxide - oxygen neighboring pair associative
        desorption; i.e., Langmuir-Hinshelwoodd reaction.

        :param self.o_co_ads: An integer that represents the minimum length of
        the state for carbon monoxide adsorption.

        :param self.o_co_des: An integer that represents the minimum length of
        the state for carbon monoxide desorption.

        :param self.o_co_dif: An integer that epresents the minimum length of the state for
        carbon monoxide diffusion.

        :param self.o_o_ads: An integer that represents the minimum length of
        the state for oxygen adsoprtion.

        :param self.o_o_des: An integer tnat represents the minimum length of
        the state for oxygen desoprtion.

        :param self.o_o_dif: An integer that represents the minimum length of
        the state for oxygen diffusion.

        Inherited parameters:

        :param self.sites: The maximum number of sites the system has.

        :param  self.states: The UNIQUE states in which each side can be in.
    """

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Getters, Setters and Deleters.
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    @property
    def k_o_ads(self):
        """ Gets the oxygen adsorption order.
        """
        return cp.deepcopy(self.__k_o_ads)

    @k_o_ads.setter
    def k_o_ads(self, _):
        """ Sets the oxygen adsorption parameter.
        """
        self.__k_o_ads = "k.O.ads"

    @k_o_ads.deleter
    def k_o_ads(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def k_o_des(self):
        """ Gets the oxygen desorption order.
        """
        return cp.deepcopy(self.__k_o_des)

    @k_o_des.setter
    def k_o_des(self, _):
        """ Sets the oxygen desorption parameter.
        """
        self.__k_o_des = "k.O.des"

    @k_o_des.deleter
    def k_o_des(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def k_o_dif(self):
        """ Gets the oxygen diffusion order.
        """
        return cp.deepcopy(self.__k_o_dif)

    @k_o_dif.setter
    def k_o_dif(self, _):
        """ Sets the oxygen diffusion parameter.
        """
        self.__k_o_dif = "k.O.dif"

    @k_o_dif.deleter
    def k_o_dif(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def k_co_ads(self):
        """ Gets the carbon monoxide adsorption order.
        """
        return cp.deepcopy(self.__k_co_ads)

    @k_co_ads.setter
    def k_co_ads(self, _):
        """ Sets the carbon monoxide adsorption parameter.
        """
        self.__k_co_ads = "k.CO.ads"

    @k_co_ads.deleter
    def k_co_ads(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def k_co_des(self):
        """ Gets the carbon monoxide desorption order.
        """
        return cp.deepcopy(self.__k_co_des)

    @k_co_des.setter
    def k_co_des(self, _):
        """ Sets the carbon monoxide desorption parameter.
        """
        self.__k_co_des = "k.CO.des"

    @k_co_des.deleter
    def k_co_des(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def k_co_dif(self):
        """ Gets the carbon monoxide diffusion order.
        """
        return cp.deepcopy(self.__k_co_dif)

    @k_co_dif.setter
    def k_co_dif(self, _):
        """ Sets the carbon monoxide diffusion parameter.
        """
        self.__k_co_dif = "k.CO.dif"

    @k_co_dif.deleter
    def k_co_dif(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def k_coo_lh(self):
        """ Gets the carbon monoxide - oxygen reaction order.
        """
        return cp.deepcopy(self.__k_coo_lh)

    @k_coo_lh.setter
    def k_coo_lh(self, _):
        """ Sets the carbon monoxide - oxygen reaction parameter.
        """
        self.__k_coo_lh = "k.COO.lh"

    @k_coo_lh.deleter
    def k_coo_lh(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def k_coo_er(self):
        """ Gets the carbon monoxide - oxygen  gas-phase reaction order.
        """
        return cp.deepcopy(self.__k_coo_er)

    @k_coo_er.setter
    def k_coo_er(self, _):
        """ Sets the carbon monoxide - oxygen  gas-phase reaction parameter.
        """
        self.__k_coo_er = "k.COO.el"

    @k_coo_er.deleter
    def k_coo_er(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_o_ads(self):
        """ Gets the oxygen adsorption order.
        """
        return cp.deepcopy(self.__o_o_ads)

    @o_o_ads.setter
    def o_o_ads(self, _):
        """ Sets the oxygen adsorption parameter.
        """
        self.__o_o_ads = 2

    @o_o_ads.deleter
    def o_o_ads(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_o_des(self):
        """ Gets the oxygen desorption order.
        """
        return cp.deepcopy(self.__o_o_des)

    @o_o_des.setter
    def o_o_des(self, _):
        """ Sets the oxygen desorption parameter.
        """
        self.__o_o_des = 2

    @o_o_des.deleter
    def o_o_des(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_o_dif(self):
        """ Gets the oxygen diffusion order.
        """
        return cp.deepcopy(self.__o_o_dif)

    @o_o_dif.setter
    def o_o_dif(self, _):
        """ Sets the oxygen diffusion parameter.
        """
        self.__o_o_dif = 2

    @o_o_dif.deleter
    def o_o_dif(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_co_ads(self):
        """ Gets the carbon monoxide adsorption order.
        """
        return cp.deepcopy(self.__o_co_ads)

    @o_co_ads.setter
    def o_co_ads(self, _):
        """ Sets the carbon monoxide adsorption parameter.
        """
        self.__o_co_ads = 1

    @o_co_ads.deleter
    def o_co_ads(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_co_des(self):
        """ Gets the carbon monoxide desorption order.
        """
        return cp.deepcopy(self.__o_co_des)

    @o_co_des.setter
    def o_co_des(self, _):
        """ Sets the carbon monoxide desorption parameter.
        """
        self.__o_co_des = 1

    @o_co_des.deleter
    def o_co_des(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_co_dif(self):
        """ Gets the carbon monoxide diffusion order.
        """
        return cp.deepcopy(self.__o_co_dif)

    @o_co_dif.setter
    def o_co_dif(self, _):
        """ Sets the carbon monoxide diffusion parameter.
        """
        self.__o_co_dif = 2

    @o_co_dif.deleter
    def o_co_dif(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_coo_lh(self):
        """ Gets the carbon monoxide - oxygen reaction order.
        """
        return cp.deepcopy(self.__o_coo_lh)

    @o_coo_lh.setter
    def o_coo_lh(self, _):
        """ Sets the carbon monoxide - oxygen reaction parameter.
        """
        self.__o_coo_lh = 1

    @o_coo_lh.deleter
    def o_coo_lh(self):
        pass

    # --------------------------------------------------------------------------

    @property
    def o_coo_er(self):
        """ Gets the carbon monoxide - oxygen  gas-phase reaction order.
        """
        return cp.deepcopy(self.__o_coo_er)

    @o_coo_er.setter
    def o_coo_er(self, _):
        """ Sets the carbon monoxide - oxygen  gas-phase reaction parameter.
        """
        self.__o_coo_er = 1

    @o_coo_er.deleter
    def o_coo_er(self):
        pass

    # --------------------------------------------------------------------------
    # Constructor.
    # --------------------------------------------------------------------------

    def __init__(self, sites=1, states=("E",)):
        """ Builds a class that writes the equations for the carbon monoxide -
            oxygen associative reaction on ruthenium (111);
            J. Chem. Phys. 143, 204702 (2015).

            Initializes the class with the standard parameters.
        """

        # Initialize the super class.
        super(COOxidation, self).__init__(sites, ('CO', 'O', 'E'))

        # ----------------------------------------------------------------------
        # Define order of the process; i.e., minimum length of the state for the
        # process to happen.
        # ----------------------------------------------------------------------

        # Represents the minimum length of the state for oxygen related
        # processes.
        self.o_o_ads = 2
        self.o_o_des = 2
        self.o_o_dif = 2

        # Represents the minimum length of the state for carbon monoxide related
        # processes.
        self.o_co_ads = 1
        self.o_co_des = 1
        self.o_co_dif = 2

        # Represents the minimum length of the state for carbon monoxide -
        # oxygen reaction related processes.
        self.o_coo_lh = 2
        self.o_coo_er = 1

        # ----------------------------------------------------------------------
        # Define the strings for the constant names.
        # ----------------------------------------------------------------------

        # String that represents the rate constants of oxygen related processes.
        self.k_o_ads = "k.O.ads"
        self.k_o_des = "k.O.des",
        self.k_o_dif = "k.O.Dif"

        # String that represents the rate constants of carbon monoxide related
        # processes.
        self.k_co_ads = "k.CO.ads"
        self.k_co_des = "k.CO.des"
        self.k_co_dif = "k.CO.dif"

        # String that represents the rate constants of carbon monoxide - oxygen
        # reaction related processes.
        self.k_coo_lh = "k.COO.lh"
        self.k_coo_er = "k.COO.er"

