# simple likelihood for sdss_lrgDR4
from likelihood_class import likelihood
import numpy as np
import math

class sdss_lrgDR4(likelihood):

  def __init__(self,path,data,command_line,log_flag,default):

    likelihood.__init__(self,path,data,command_line,log_flag,default)

    # require P(k) from class
    self.need_cosmo_arguments(data,{'output':'mPk'})
 
    # In case of a comparison, stop here
    if not default:
      return

    # in this likelihood, the P(k,z) is supposed to be evaluated only in z=0 
    self.z_size=1
    self.z = np.zeros((self.z_size),'float64')
    self.z[0]=0

    # read values of k (in h/Mpc)
    self.k_size=self.max_mpk_kbands_use-self.min_mpk_kbands_use+1
    self.mu_size=1
    self.k = np.zeros((self.k_size,self.z_size,self.mu_size),'float64')
    self.kh = np.zeros((self.k_size,self.z_size,self.mu_size),'float64')

    datafile = open(self.data_directory+self.kbands_file,'r')
    for i in range(self.num_mpk_kbands_full):
      line = datafile.readline()
      if i+2 > self.min_mpk_kbands_use and i < self.max_mpk_kbands_use:
        self.kh[i-self.min_mpk_kbands_use+1,0,0]=float(line.split()[0])
    datafile.close()      

    khmax = self.kh[-1,0,0]
    
    # require k_max from the cosmological module
    self.need_cosmo_arguments(data,{'P_k_max_h/Mpc':khmax})
 
    # In case of a comparison, stop here
    if not default:
      return

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

    # reduced Hubble parameter
    h = _cosmo._h()
    
    # redhsift used in the rescaling
    z_rescaling = 0.35
 
    # angular diameter distance at this redshift, in Mpc/h
    d_angular = _cosmo._angular_distance(z_rescaling)
    d_angular *= h  

    # radial distance at this redshift, in Mpc/h
    r,Hz = _cosmo.z_of_r([z_rescaling])
    d_radial = z_rescaling*h/Hz[0]

    # scaling factor = (d_angular**2 * d_radial)^(1/3) relative
    # to a flat reference model with Omega_m=0.25, Omega_Lambda=0.75
    scaling = pow((d_angular/722.4134)**2*(d_radial/897.9993),1./3.)
    #scaling = pow((d_angular/1032.0191)**2*(d_radial/1282.8561),1./3.)

    # get rescaled values of k in 1/Mpc
    self.k=self.kh *h*scaling

    # get P(k) at right values of k, convert it to (Mpc/h)^3 and rescale it
    P_lin = np.zeros((self.k_size,self.z_size,self.mu_size),'float64')
    P_lin = _cosmo._get_pk(self.k,self.z,self.k_size,self.z_size,self.mu_size)
    P_lin /= (h*scaling)**3

    do_marge = self.Q_marge

    if do_marge and self.Q_flat:

      P_th =  np.zeros((self.k_size),'float64')
      for i in range(self.k_size):
        P_th[i] = P_lin[i,0,0]/(1.+self.Ag*self.k[i]) 

      k2 =  np.zeros((self.k_size),'float64')
      for i in range(self.k_size):
        k2[i] = P_th[i] * self.k[i]**2

      W_P_th =  np.zeros((self.n_size),'float64')
      W_P_th = np.dot(self.window,P_th)

      W_P_th_k2 =  np.zeros((self.n_size),'float64')
      W_P_th_k2 = np.dot(self.window,k2)
                
      w = np.zeros((self.n_size),'float64') 
      w=1./(self.P_err*self.P_err)

      covdat = np.zeros((self.n_size),'float64') 
      covdat=self.P_obs*w  

      covth = np.zeros((self.n_size),'float64') 
      covth=W_P_th*w

      covth_k2 = np.zeros((self.n_size),'float64') 
      covth_k2=W_P_th_k2*w

      offdiag=sum(covth*W_P_th_k2)

      Mat=np.zeros((2,2),'float64')
      Mat=[[sum(covth*W_P_th),offdiag],[offdiag,sum(covth_k2*W_P_th_k2)]]

      Vec=np.zeros((2),'float64')
      Vec=[sum(covdat*W_P_th),sum(covdat*W_P_th_k2)]

      chi2=-sum(self.P_obs*covdat)+np.dot(Vec,np.dot(np.linalg.inv(Mat),Vec))-math.log(np.linalg.det(Mat))

    else:  

      print "case without marginalization and Q_flat not yet coded"
      exit()

    return - chi2/2.
