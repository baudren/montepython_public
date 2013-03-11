import os
import numpy as np
from likelihood_class import likelihood

class da_rec(likelihood):
  
  # The initialization routine is no different from the parent class one.

  # compute likelihood
  def loglkl(self, cosmo, data):

    da = cosmo._angular_distance(self.z_rec)
    chi2 = ((da-self.da_rec)/self.da_error)**2

    # return ln(L)
    lkl = -0.5*chi2 

    return lkl
