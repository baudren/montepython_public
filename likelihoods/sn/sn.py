import os
import numpy as np
from likelihood_class import likelihood

class sn(likelihood):
  
  def __init__(self,path,command_line=False):

    likelihood.__init__(self,path,command_line)

    #abs_folder = os.path.dirname(os.path.abspath(path))+'/'
    

    self.z	= np.array([],'float64')
    self.moduli = np.array([],'float64')
    for line in open(self.z_mu_dmu,'r'):
      self.z		= np.append(self.z,float(line.split()[1]))
      self.moduli 	= np.append(self.z,float(line.split()[2]))


    if self.has_syscovmat:
      covmat_filename = self.covmat_sys
    else:
      covmat_filename = self.covmat_nosys

    covmat = np.zeros((np.shape(self.z)[0],np.shape(self.z)[0]),'float64')
    i=0
    for line in open(covmat_filename,'r'):
      covmat[i] = line.split()
      i+=1

    self.inv_covmat=np.linalg.inv(covmat)
    self.trace_inv_covmat=np.trace(self.inv_covmat)

  def _loglkl(self,_cosmo,data):

    difference = np.ndarray(np.shape(self.z)[0],'float64')

    for i in range(np.shape(self.z)[0]):
      d = _cosmo._angular_distance(self.z[i])
      #d = 2*self.z[i]
      difference[i] = 5.0*np.log((1.+ self.z[i] )**2*d)/np.log(10) + 25.0 - self.moduli[i]

    AT = np.dot(np.transpose(difference),np.dot(self.inv_covmat,difference))
    if self.has_marginalization:
      BT = np.sum(np.dot(self.inv_covmat,difference))
    else:
      BT = 0

    chi_squared = AT - (BT**2)/self.trace_inv_covmat
    self.loglkl = - 0.5 * chi_squared 
    return self.loglkl
