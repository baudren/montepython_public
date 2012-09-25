# simple likelihood for sdss_lrgDR4
from likelihood_class import likelihood
import numpy as np

class sdss_lrgDR4(likelihood):

  def __init__(self,path,data,command_line,log_flag,default):

    likelihood.__init__(self,path,data,command_line,log_flag,default)

    self.need_Class_arguments(data,{'output':'mPk', 'P_k_max_h/Mpc':2})
 
    # In case of a comparison, stop here
    if not default:
      return

    # read values of k (in h/Mpc) and z
    self.k_size=self.max_mpk_kbands_use-self.min_mpk_kbands_use+1
    self.z_size=1
    self.mu_size=1

    self.z = np.zeros((self.z_size),'float64')
    self.z[0]=0

    self.k = np.zeros((self.k_size,self.z_size,self.mu_size),'float64')
    self.kh = np.zeros((self.k_size,self.z_size,self.mu_size),'float64')

    datafile = open(self.data_directory+self.kbands_file,'r')
    for i in range(self.num_mpk_kbands_full):
      line = datafile.readline()
      if i+2 > self.min_mpk_kbands_use and i < self.max_mpk_kbands_use:
          self.kh[i-self.min_mpk_kbands_use+1,0,0]=float(line.split()[0])
    datafile.close()      
    
    # read window functions
    self.n_size=self.max_mpk_points_use-self.min_mpk_points_use+1

    self.window=np.zeros((self.n_size,self.k_size),'float64')

    datafile = open(self.data_directory+self.windows_file,'r')
    for i in range(self.num_mpk_points_full):
      line = datafile.readline()
      if i+2 > self.min_mpk_points_use and i < self.max_mpk_points_use:
        for j in range(self.k_size):  
          self.window[i-self.min_mpk_points_use+1,j]=float(line.split()[j+self.min_mpk_kbands_use-1])
    datafile.close()    

    # read measurements
    self.P_obs=np.zeros((self.n_size),'float64')
    self.P_err=np.zeros((self.n_size),'float64')

    datafile = open(self.data_directory+self.measurements_file,'r')
    for i in range(2):
      line = datafile.readline()
    for i in range(self.num_mpk_points_full):
      line = datafile.readline()
      if i+2 > self.min_mpk_points_use and i < self.max_mpk_points_use:
        self.P_obs[i-self.min_mpk_points_use+1]=float(line.split()[3])  
        self.P_err[i-self.min_mpk_points_use+1]=float(line.split()[4])  
    datafile.close()

    return

  # compute likelihood

  def loglkl(self,_cosmo,data):

    # get values of k in 1/Mpc
    self.k=self.kh *_cosmo._h()

    # get P(k) at right values of k and convert it to (Mpc/h)^3
    P_th = np.zeros((self.k_size,self.z_size,self.mu_size),'float64')
    P_th = _cosmo._get_pk(self.k,self.z,self.k_size,self.z_size,self.mu_size)
    P_th /= (_cosmo._h())**3

    # intermediate quantitites for chi2    
    W_P_th =  np.zeros((self.n_size),'float64')
    for i in range(self.n_size):
      W_P_th[i]=0
      for j in range(self.k_size):
        W_P_th += self.window[i,j] * P_th[j,0,0]

    w = np.zeros((self.n_size),'float64') 
    normV=0.
    tmp=0.
    for i in range(self.n_size):   
      w[i] = 1/(self.P_err[i]*self.P_err[i])
      normV += W_P_th[i]*W_P_th[i]*w[i]
      tmp += W_P_th[i]*self.P_obs[i]*w[i]
    tmp /= normV  

    # final chi2 expression (marginalized over bias)
    chi2 = 0.
    for i in range(self.n_size):
      chi2 += self.P_obs[i]*(self.P_obs[i]-W_P_th[i]*tmp)*w[i]  

    #print 'sdss chi2=',chi2  

    return - chi2/2.
