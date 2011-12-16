import os
import numpy as np
from likelihood_class import likelihood

class sn(likelihood):
  
  # initialization routine

  def __init__(self,path,data,command_line=False):

    likelihood.__init__(self,path,data,command_line)

    # define array for values of z and data points
    self.z	= np.array([],'float64')
    self.moduli = np.array([],'float64')

    # read redshifts and data points
    for line in open(self.z_mu_dmu,'r'):
      self.z		= np.append(self.z,float(line.split()[1]))
      self.moduli 	= np.append(self.moduli,float(line.split()[2]))

    # number of data points
    self.num_points=np.shape(self.z)[0]

    # define correlation m,atrix
    covmat = np.zeros((self.num_points,self.num_points),'float64')

    # file containing correlation matrix
    if self.has_syscovmat:
      covmat_filename = self.covmat_sys
    else:
      covmat_filename = self.covmat_nosys

    # read correlation matrix
    i=0
    for line in open(covmat_filename,'r'):
      covmat[i] = line.split()
      i+=1

    # invert correlation matrix
    self.inv_covmat=np.linalg.inv(covmat)

    # find sum of all matrix elements (sounds odd that there is 
    # not a trace here instead, but this is correct!)
    self.inv_covmat_sum=np.sum(self.inv_covmat)

    # end of initialization

  # compute likelihood

  def _loglkl(self,_cosmo,data):

    # define array for difference between theory and observations
    difference = np.ndarray(self.num_points,'float64')

    # for each point, compute luminosity distance d_L=(1+z)**2d_A and infer
    # theoretical prediction and difference with observation
    for i in range(self.num_points):
      d = _cosmo._angular_distance(self.z[i])
      difference[i] = 5*np.log((1+self.z[i])**2*d)/np.log(10) + 25 - self.moduli[i]

    # chisquare before analytic marginalization
    AT = np.dot(difference,np.dot(self.inv_covmat,difference))

    # correct chi square with effect of analytic marginalization
    if self.has_marginalization:
      BT = np.sum(np.dot(self.inv_covmat,difference))
    else:
      BT = 0

    # final chi square
    chi_squared = AT - (BT**2)/self.inv_covmat_sum

    # return ln(L)
    self.loglkl = - 0.5 * chi_squared 

    return self.loglkl
