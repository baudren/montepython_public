# simple likelihood for WiggleZ
from likelihood_class import likelihood
import numpy as np
import math
import os

class WiggleZ(likelihood):

  def __init__(self,path,data,command_line,log_flag,default):

    likelihood.__init__(self,path,data,command_line,log_flag,default)

    # require P(k) from class
    self.need_cosmo_arguments(data,{'output':'mPk'})
 
    # in this likelihood, the P(k,z) is evaluated in 4 redshift bins, 
    # 0.1 < z < 0.3, 0.3 < z < 0.5, 0.5 < z < 0.7, 0.7 < z < 0.9
    self.z_size  = 4
    self.z_mean  = np.zeros((self.z_size),  'float64')
    self.z_edges = np.zeros((self.z_size+1),'float64')
    self.z_mean  = [0.2,0.4,0.6,0.8]
    self.z_eges  = [0.1,0.3,0.5,0.7,0.9]

    # read all values of k (in h/Mpc). Note that for a given redshift bin, not
    # all of these values will be used
    self.k = np.zeros((self.num_mpk_kbands_full),'float64')
    self.kh = np.zeros((self.num_mpk_kbands_full),'float64')

    datafile = open(self.data_directory+self.kbands_file,'r')
    for i in range(self.num_mpk_kbands_full):
      line = datafile.readline()
      self.kh[i]=float(line.split()[0])
    datafile.close()      

    khmax = self.kh[-1]
    
    # require k_max from the cosmological module
    self.need_cosmo_arguments(data,{'P_k_max_h/Mpc':khmax})
 
    # In case of a comparison, stop here
    if not default:
      return

    # read all .dataset files, and assing it to a particular dictionnary
    self.redshift_bins = []
    for elem in self.dataset_files:
      self.read_dataset_file(elem)
    exit()

    return

  def read_dataset_file(self,File):
    # Opening the dataset file, that contains, for each redshift bin, all relevant quantities.
    print 'reading',File

    Dict = {}

    # Reading the dataset file
    if os.path.isfile(self.data_directory+File):
      dataset = open(self.data_directory+File,'r')
      for line in dataset:
        if (line.find('#') == -1 and len(line)>3):
          result = str(line.split('=')[-1].replace(' ','').replace('\n',''))
          # Eliminating the data/ part of the filenames
          if result.find('data/') != -1:
            result = result.replace('data/','')
          # Try converting to a float or an int. If no success, then it really is a string
          if result.find('.') != -1:
            try: result = float(result)
            except ValueError: pass
          else:
            try: result = int(result)
            except ValueError: pass
          if result == 'F':
            result = False
          elif result == 'T':
            result = True
          # Recovering the name of the field
          name   = str(line.split('=')[0].replace(' ',''))
          Dict[name] = result

    dataset.seek(0)
    dataset.close()
    
    # read window functions
    Dict['n_size'] = Dict['max_mpk_points_use']-Dict['min_mpk_points_use']+1
    Dict['k_size'] = Dict['max_mpk_kbands_use']-Dict['min_mpk_kbands_use']+1

    Dict['window'] = np.zeros((Dict['n_size'],Dict['k_size']),'float64')
    datafile = open(self.data_directory+Dict['windows_file'],'r')
    for i in range(Dict['num_mpk_kbands_full']):
      line = datafile.readline()
      if i+2 > Dict['min_mpk_points_use'] and i < Dict['max_mpk_points_use']:
        for j in range(Dict['k_size']):  
          Dict['window'][i-Dict['min_mpk_points_use']+1,j]=float(line.split()[j+Dict['min_mpk_kbands_use']-1])
    datafile.close()    

    # read measurements
    Dict['P_obs'] = np.zeros((Dict['n_size']),'float64')
    Dict['P_err'] = np.zeros((Dict['n_size']),'float64')

    datafile = open(self.data_directory+Dict['measurements_file'],'r')
    # WARNING: This is not robust if presentation changes
    for i in range(2):
      line = datafile.readline()
    for i in range(Dict['num_mpk_points_full']):
      line = datafile.readline()
      if i+2 > Dict['min_mpk_points_use'] and i < Dict['max_mpk_points_use']:
        Dict['P_obs'][i-Dict['min_mpk_points_use']+1]=float(line.split()[3])  
        Dict['P_err'][i-Dict['min_mpk_points_use']+1]=float(line.split()[4])  
    datafile.close()

    # Finally, appending the dictionary of values for this redshift bin to the list redshift_bins
    self.redshift_bins.append(Dict)
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
