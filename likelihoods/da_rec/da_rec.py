import os
import numpy as np
from likelihood_class import likelihood

class da_rec(likelihood):
  
  # initialization routine

  def __init__(self,path,data,command_line,log_flag,default):

    likelihood.__init__(self,path,data,command_line,log_flag,default)
    if not default:
      return
    # end of initialization

  # compute likelihood

  def loglkl(self,_cosmo,data):

    da = _cosmo._angular_distance(self.z_rec)
    chi2 = ((da-self.da_rec)/self.da_error)**2

    # return ln(L)
    lkl = - 0.5 * chi2 

    return lkl
