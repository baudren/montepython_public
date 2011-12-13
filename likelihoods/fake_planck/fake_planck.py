import os
import numpy as np
try:
  import clik
except ImportError:
  print " /|\  You must first activate the binaries from the Clik distribution,"
  print "/_o_\ please run : source /path/to/clik/bin/clik_profile.sh"
  print "      and try again."
  exit()
from likelihood_class import likelihood

class fake_planck(likelihood):

  def __init__(self,path,command_line=False):
    likelihood.__init__(self,path,command_line)
    self.clik = clik.clik(self.path_clik)

  def _loglkl(self,_cosmo,data):

    T = _cosmo._T_cmb()
    T_squared=np.array(T**2.*1e12,'float64')
    cl=np.array(_cosmo.lensed_cl(),'float64')
    tot=np.zeros(np.sum(self.clik.lmax)+6)
    for i in range(cl.shape[0]-1):
      ieq=i
      for j in range(1701):
	bound=1701
	tot[j+i*bound]=np.array(cl[ieq][j]*T_squared,'float64')
    loglkl=self.clik(tot)[0]
    return loglkl
