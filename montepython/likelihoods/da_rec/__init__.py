import os
import numpy as np
from montepython.likelihood_class import Likelihood


class da_rec(Likelihood):

    # The initialization routine is no different from the parent class one.

    # compute likelihood
    def loglkl(self, cosmo, data):

        da = cosmo.angular_distance(self.z_rec)
        chi2 = ((da - self.da_rec) / self.da_error) ** 2

        # return ln(L)
        lkl = -0.5 * chi2

        return lkl
