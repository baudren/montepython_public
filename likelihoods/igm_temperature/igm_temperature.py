import os
import numpy as np
from likelihood_class import likelihood_prior

class igm_temperature(likelihood_prior):
  
  def loglkl(self, cosmo, data):

    lkl = 0

    if(cosmo._baryon_temperature(self.z)>self.Tb_max):
      lkl=data.boundary_loglike
      #print 'Tb(',self.z,') =',self.cosmo._igm_temperature(self.z),' > ',self.Tb_max

    return lkl
