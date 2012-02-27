from likelihood_class import likelihood
import os
import numpy as np
import math
# Adapted from JL

class euclid_pk(likelihood):

  def __init__(self,path,data,command_line,log_flag,default):

    likelihood.__init__(self,path,data,command_line,log_flag,default)

    self.need_Class_arguments(data,{'output':'mPk'})
    
    #################
    # find number of galaxies for each mean redshift value
    #################

    # Compute the number of galaxies for each \bar z
    # For this, one needs dn/dz TODO
    # then n_g(\bar z) = int_{\bar z - dz/2}^{\bar z + dz/2} dn/dz dz
    # self.z_mean will contain the central values

    self.n_g = np.zeros(self.nbin,'float64')
    self.z_mean = np.zeros(self.nbin,'float64')
    i=0
    for z in np.arange(self.zmin,self.zmax+self.dz,self.dz):
      self.z_mean[i] = z
      i+=1

    # For each bin, compute the biais function,

    self.b = np.zeros(self.nbin,'float64')
    for Bin in range(self.nbin):
      self.b[Bin] = math.sqrt(self.z_mean[Bin]+1.)

    self.need_Class_arguments(data,{'z_max_pk':self.z_mean[-1]+self.dz})
    # Force Class to store Pk for k up to an arbitrary number (since self.r is not yet decided)... TODO
    self.need_Class_arguments(data,{'P_k_max_1/Mpc':5})

    # In case of a comparison, stop here
    if not default:
      return

    ################
    # Noise spectrum  TODO
    ################

    # If the file exists, initialize the fiducial values
    #self.Cl_fid = np.zeros((self.nlmax,self.nbin,self.nbin),'float64')
    self.fid_values_exists = False
    self.H_fid   = np.zeros(self.nbin,'float64')
    self.D_A_fid = np.zeros(self.nbin,'float64') 
    if os.path.exists(self.data_directory+'/'+self.fiducial_file):
      self.fid_values_exists = True
      fid_file = open(self.data_directory+'/'+self.fiducial_file,'r')
      line = fid_file.readline()
      while line.find('#')!=-1:
	line = fid_file.readline()
      while (line.find('\n')!=-1 and len(line)==1):
	line = fid_file.readline()
      for nl in range(self.nlmax):
	for Bin in range(self.nbin):
	  for Bin2 in range(self.nbin):
	    self.Cl_fid[nl,Bin,Bin2] = float(line)
	    line = fid_file.readline()
      
    # Else the file will be created in the loglkl() function. 
    return
  
  # Maybe to remove ? TODO
  def galaxy_distribution(self,z):

    zmean = 0.9
    z0 = zmean/1.412

    galaxy_dist = z**2*math.exp(-(z/z0)**(1.5))
    
    return galaxy_dist

  def loglkl(self,_cosmo,data):
    # First thing, recover the angular distance and Hubble factor for each
    # redshift
    H   = np.zeros(self.nbin,'float64')
    D_A = np.zeros(self.nbin,'float64') 
    r   = np.zeros(self.nbin,'float64')

    # H is incidentally also dz/dr
    r,H = _cosmo.z_of_r(self.z_mean)
    D_A = _cosmo._angular_distance(self.z_mean)
    #for Bin in range(self.nbin):
      #print '%.4g %.4g %.4g' % (r[Bin],H[Bin],D_A[Bin])

    # Compute sigma_r = dr(z)/dz sigma_z with sigma_z = 0.001(1+z)
    sigma_r = np.zeros(self.nbin,'float64')
    for Bin in range(self.nbin):
      sigma_r[Bin] = 0.001*(1.+self.z_mean[Bin])/H[Bin]

    # Compute V_survey, for each given redshift bin, which is the volume of a
    # shell times the sky coverage:
    V_survey = np.zeros(self.nbin,'float64')
    for Bin in range(self.nbin):
      V_survey[Bin] = 4.*math.pi*(r[Bin]**2)*(1+self.z_mean[Bin])**-3 *self.dz/H[Bin]
    print V_survey
    return 0
