import os
import numpy as np
from montepython.likelihood_class import Likelihood_prior


class igm_temperature(Likelihood_prior):

    def loglkl(self, cosmo, data):

        lkl = 0

        if(cosmo.baryon_temperature(self.z) > self.Tb_max):
            lkl = data.boundary_loglike
            # print 'Tb(',self.z,') =',self.cosmo.igm_temperature(self.z),' >
            # ',self.Tb_max

        return lkl
