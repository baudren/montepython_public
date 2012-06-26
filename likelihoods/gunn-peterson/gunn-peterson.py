import os
import numpy as np
from likelihood_class import likelihood

class gunn_peterson(likelihood_prior):
  
  def loglkl(self,_cosmo,data):

    xHI_reio = 1.-_cosmo._ionization_fraction(self.z_reio)
    xHI_noreio = 1.-_cosmo._ionization_fraction(self.z_noreio)

    lkl=0

    if (xHI_reio>gunn_peterson.xHI_max) lkl=1.e10
    if (xHI_noreio<gunn_peterson.xHI_min) lkl=1.e10

    print "at z=5.5:"xHI_reio,", at z=6:",xHI_reio,lkl

    return lkl
