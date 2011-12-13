import os
import numpy as np
try:
  import clik
except ImportError:
  print " /|\  You must first activate the binaries from the Clik distribution,"
  print "/_o_\ please run : source /path/to/clik/bin/clik_profile.sh"
  print "      and try again."
  exit()

class fake_planck():

  def __init__(self,path='fake_planck.data'):

    abs_folder = os.path.dirname(os.path.abspath(path))+'/'

    for line in open(path,'r'):
      if line.find('#')==-1:
	if line.find('fake_planck.')!=-1:
	  exec(line.replace('fake_planck.','self.'))
    
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
