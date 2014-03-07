import os
import numpy as np
from montepython.likelihood_class import Likelihood_prior


class gunn_peterson(Likelihood_prior):

    def loglkl(self, cosmo, data):

        xHI_reio = 1. - cosmo.ionization_fraction(self.z_reio)
        xHI_noreio = 1. - cosmo.ionization_fraction(self.z_noreio)

        lkl = 0

        if (xHI_reio > self.xHI_max):
            lkl = data.boundary_loglike
        if (xHI_noreio < self.xHI_min):
            lkl = data.boundary_loglike

        return lkl
