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

class clik_wmap_full(likelihood):

  def __init__(self,path,data,command_line=False):
    likelihood.__init__(self,path,data,command_line)
    self._need_Class_args(data,{'lensing':'yes', 'output':'tCl lCl pCl'})
    self.clik = clik.clik(self.path_clik)

  def _loglkl(self,_cosmo,data):

    T = _cosmo._T_cmb()
    T_squared=np.array(T**2.*1e12,'float64')
    cl=np.array(_cosmo.lensed_cl(),'float64')
    tot=np.zeros(np.sum(self.clik.get_lmax())+6)
    for i in range(cl.shape[0]):
      ieq=i
      for j in range(1201):
	bound=1201
	tot[j+i*bound]=np.array(cl[ieq][j]*T_squared,'float64')
    loglkl=self.clik(tot)[0]
    return loglkl

