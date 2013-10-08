from likelihood_class import likelihood
import os
import numpy as np
from math import exp,log,sqrt,pi

class euclid_pk(likelihood):

  def __init__(self, path, data, command_line):

    likelihood.__init__(self, path, data, command_line)

    self.need_cosmo_arguments(data, {'output':'mPk'})
    
    #################
    # find number of galaxies for each mean redshift value
    #################

    # Compute the number of galaxies for each \bar z
    # For this, one needs dn/dz TODO
    # then n_g(\bar z) = int_{\bar z - dz/2}^{\bar z + dz/2} dn/dz dz
    # self.z_mean will contain the central values

    self.z_mean = np.zeros(self.nbin,'float64')

    # Deduce the dz step from the number of bins and the edge values of z
    self.dz = (self.zmax-self.zmin)/(self.nbin-1.)
    i=0
    for z in np.arange(self.zmin,self.zmax+self.dz,self.dz):
      self.z_mean[i] = z
      i+=1

    # Store the z edge values
    self.z_edges = np.zeros(self.nbin+1,'float64')
    i = 0
    for z in np.arange(self.zmin-self.dz/2.,self.zmax+self.dz,self.dz):
      self.z_edges[i] = z
      i+=1

    # Store the total vector z, with edges + mean
    self.z = np.zeros(2*self.nbin+1,'float64')
    i=0
    for z in np.arange(self.zmin-self.dz/2.,self.zmax+self.dz,self.dz/2.):
      self.z[i] = z
      i+=1

    self.need_cosmo_arguments(data,{'z_max_pk':self.z[-1]})

    # For each bin, compute the biais function,
    self.b = np.zeros(self.nbin,'float64')
    for Bin in range(self.nbin):
      self.b[Bin] = sqrt(self.z_mean[Bin]+1.)

    # Force the cosmological module to store Pk for k up to an arbitrary number
    # (since self.r is not yet decided)... TODO
    self.need_cosmo_arguments(data,{'P_k_max_1/Mpc':1.5*self.kmax})

    # Define the k values for the integration (from kmin to kmax), at which the
    # spectrum will be computed (and stored for the fiducial model)
    # k_size is deeply arbitrary here, TODO
    
    self.k_fid = np.zeros(self.k_size,'float64')
    for i in range(self.k_size):
      self.k_fid[i] = exp( i*1.0 /(self.k_size-1) * log(self.kmax/self.kmin) + log(self.kmin))

    ################
    # Noise spectrum TODO properly
    ################

    self.n_g = np.zeros(self.nbin,'float64')
    #for index_z in range(self.nbin):
      #print self.z_mean[index_z],self.z[2*index_z],self.z[2*index_z+2]
      #self.n_g[index_z] = self.galaxy_distribution(self.z[2*index_z+2]) - self.galaxy_distribution(self.z[2*index_z])
    #print self.n_g

    self.n_g[0] = 6844.945
    self.n_g[1] = 7129.45
    self.n_g[2] = 7249.912
    self.n_g[3] = 7261.722
    self.n_g[4] = 7203.825
    self.n_g[5] = 7103.047
    self.n_g[6] = 6977.571
    self.n_g[7] = 6839.546
    self.n_g[8] = 6696.957
    self.n_g[9] = 5496.988
    self.n_g[10] = 4459.240
    self.n_g[11] = 3577.143
    self.n_g[12] = 2838.767
    self.n_g[13] = 2229.282
    self.n_g[14] = 1732.706
    self.n_g[15] = 1333.091

    self.n_g = self.n_g * self.efficiency * 41253.

    # If the file exists, initialize the fiducial values,
    # the spectrum will be read first, with k_size values of k and nbin values
    # of z. Then, H_fid and D_A fid will be read (each with nbin values).
    #self.Cl_fid = np.zeros((self.nlmax,self.nbin,self.nbin),'float64')
    self.fid_values_exist = False
    self.pk_nl_fid = np.zeros((self.k_size,2*self.nbin+1),'float64')
    self.H_fid       = np.zeros(2*self.nbin+1,'float64')
    self.D_A_fid     = np.zeros(2*self.nbin+1,'float64') 
    self.sigma_r_fid = np.zeros(self.nbin,'float64')

    if os.path.exists(self.data_directory+'/'+self.fiducial_file):
      self.fid_values_exist = True
      fid_file = open(self.data_directory+'/'+self.fiducial_file,'r')
      line = fid_file.readline()
      while line.find('#')!=-1:
	line = fid_file.readline()
      while (line.find('\n')!=-1 and len(line)==1):
	line = fid_file.readline()
      for index_k in range(self.k_size):
	for index_z in range(2*self.nbin+1):
	  self.pk_nl_fid[index_k,index_z] = float(line)
	  line = fid_file.readline()
      for index_z in range(2*self.nbin+1):
	self.H_fid[index_z]   = float(line.split()[0])
	self.D_A_fid[index_z] = float(line.split()[1])
	line = fid_file.readline()
      for index_z in range(self.nbin):
	self.sigma_r_fid[index_z] = float(line)
	line = fid_file.readline()
      fid_file.seek(0)
      fid_file.close()
      
    # Else the file will be created in the loglkl() function. 
    return
  
  # Galaxy distribution, returns the function D(z) from the notes
  def galaxy_distribution(self,z):

    zmean = 0.9
    z0 = zmean/1.412

    galaxy_dist = z**2*exp(-(z/z0)**(1.5))
    
    return galaxy_dist

  def loglkl(self, cosmo, data):
    # First thing, recover the angular distance and Hubble factor for each
    # redshift
    H   = np.zeros(2*self.nbin+1,'float64')
    D_A = np.zeros(2*self.nbin+1,'float64') 
    r   = np.zeros(2*self.nbin+1,'float64')

    # H is incidentally also dz/dr
    r,H = cosmo.z_of_r(self.z)
    for i in range(len(D_A)):
      D_A[i] = cosmo.angular_distance(self.z[i])

    # Compute sigma_r = dr(z)/dz sigma_z with sigma_z = 0.001(1+z)
    sigma_r = np.zeros(self.nbin,'float64')
    for index_z in range(self.nbin):
      sigma_r[index_z] = 0.001*(1.+self.z_mean[index_z])/H[2*index_z+1]

    # Compute V_survey, for each given redshift bin, which is the volume of a
    # shell times the sky coverage:
    self.V_survey = np.zeros(self.nbin,'float64')
    for index_z in range(self.nbin):
      self.V_survey[index_z] = 4.*pi*self.fsky*(r[2*index_z+1]**2)*(1+self.z_mean[index_z])**(-3) *self.dz/H[2*index_z+1]
    
    # Define the mu scale
    mu = np.zeros(self.mu_size,'float64')
    for index_mu in range(self.mu_size):
      mu[index_mu] = 2.*(index_mu)/(self.mu_size-1) - 1.

    # If the fiducial model does not exists, recover the power spectrum and
    # store it, then exit.
    if self.fid_values_exist is False:
      pk = np.zeros((self.k_size,2*self.nbin+1),'float64')
      fid_file = open(self.data_directory+'/'+self.fiducial_file,'w')
      fid_file.write('# Fiducial parameters')
      for key,value in data.mcmc_parameters.iteritems():
	fid_file.write(', %s = %.5g' % (key,value['current']*value['scale']))
      fid_file.write('\n')
      for index_k in range(self.k_size):
	for index_z in range(2*self.nbin+1):
	  pk[index_k,index_z] = cosmo.pk(self.k_fid[index_k],self.z[index_z])
	  fid_file.write('%.8g\n' % pk[index_k,index_z])
      for index_z in range(2*self.nbin+1):
	fid_file.write('%.8g %.8g\n' % (H[index_z],D_A[index_z]))
      for index_z in range(self.nbin):
	fid_file.write('%.8g\n' % sigma_r[index_z])
      print '\n\n /|\  Writting fiducial model in {0}'.format(self.data_directory+self.fiducial_file)
      print '/_o_\ for {0} likelihood'.format(self.name)
      return 1

    # NOTE: Many following loops will be hidden in a very specific numpy
    # expression, for (a more than significant) speed-up. All the following
    # loops keep the same pattern.  The colon denotes the whole range of
    # indices, so beta_fid[:,index_z] denotes the array of length self.k_size
    # at redshift z[index_z]

    # Compute the beta_fid function, for observed spectrum,
    # beta_fid(k_fid,z) = 1/2b(z) * d log(P_nl_fid(k_fid,z))/d log a
    #                   = -1/2b(z)* (1+z) d log(P_nl_fid(k_fid,z))/dz 
    
    beta_fid = np.zeros((self.k_size,self.nbin),'float64')
    for index_z in range(self.nbin):
      beta_fid[:,index_z] = -1./(2.*self.b[index_z]) * (1.+self.z_mean[index_z]) * np.log(self.pk_nl_fid[:,2*index_z+2] / self.pk_nl_fid[:,2*index_z])/(self.dz)

    # Compute the tilde P_fid(k_ref,z,mu) = H_fid(z)/D_A_fid(z)**2 ( 1 + beta_fid(k_fid,z)mu^2)^2 P_nl_fid(k_fid,z)exp ( -k_fid^2 mu^2 sigma_r_fid^2)
    self.tilde_P_fid = np.zeros((self.k_size,self.nbin,self.mu_size),'float64')
    for index_z in range(self.nbin):
      for index_mu in range(self.mu_size):
          self.tilde_P_fid[:,index_z,index_mu] = self.H_fid[2*index_z+1]/(self.D_A_fid[2*index_z+1]**2) * (1. + beta_fid[:,index_z]*(mu[index_mu]**2))**2 * self.pk_nl_fid[:,2*index_z+1]*np.exp( - self.k_fid[:]**2 *mu[index_mu]**2*self.sigma_r_fid[index_z]**2)

    ######################
    # TH PART
    ###################### 
    # Compute values of k based on k_fid (ref in paper), with formula (33 has to be corrected):
    # k^2 = ( (1-mu^2) D_A_fid(z)^2/D_A(z)^2 + mu^2 H(z)^2/H_fid(z)^2) k_fid ^ 2
    # So k = k (k_ref,z,mu)
    self.k = np.zeros((self.k_size,2*self.nbin+1,self.mu_size),'float64')
    for index_k in range(self.k_size):
      for index_z in range(2*self.nbin+1):
        self.k[index_k,index_z,:] = np.sqrt((1.-mu[:]**2)*self.D_A_fid[index_z]**2/D_A[index_z]**2 + mu[:]**2*H[index_z]**2/self.H_fid[index_z]**2 )*self.k_fid[index_k]


    # Recover the non-linear power spectrum from the cosmological module on all
    # the z_boundaries, to compute afterwards beta. This is pk_nl_th from the
    # notes.
    pk_nl_th = np.zeros((self.k_size,2*self.nbin+1,self.mu_size),'float64')
    pk_nl_th = cosmo.get_pk(self.k,self.z,self.k_size,2*self.nbin+1,self.mu_size)

    #for index_k in range(self.k_size):
    #  print self.k[index_k,0,0],pk_nl_th[index_k,0,0]
    #exit()

    # Recover the non_linear scale computed by halofit. If no scale was
    # affected, set the scale to one, and make sure that the nuisance parameter
    # epsilon is set to zero
    k_sigma = np.zeros(2.*self.nbin+1, 'float64')
    if (cosmo.nonlinear_method == 0):
      k_sigma[:]=1.e6
    else :
      k_sigma = cosmo.nonlinear_scale(self.z,2*self.nbin+1)

    # Define the alpha function, that will characterize the theoretical
    # uncertainty. 
    self.alpha = np.zeros((self.k_size,2*self.nbin+1,self.mu_size),'float64')
    for index_z in range(2*self.nbin+1):
      for index_mu in range(self.mu_size):
        self.alpha[:,index_z,index_mu] = np.log(1. + self.k[:,index_z,index_mu]/k_sigma[index_z]) / (1. + np.log(1.+ self.k[:,index_z,index_mu]/k_sigma[index_z]))*self.theoretical_error

    # recover the e_th part of the error function
    e_th = self.coefficient_f_nu*cosmo.Omega_nu/cosmo.Omega_m

    # Compute the Error E_th function
    #E_th = np.zeros((self.k_size,2*self.nbin+1,self.mu_size),'float64')
    #for index_z in range(2*self.nbin+1):
    #  for index_mu in range(self.mu_size):
    #    E_th[:,index_z,index_mu] = np.log(1. + self.k[:,index_z,index_mu]/k_sigma[index_z]) / (1. + np.log(1.+ self.k[:,index_z,index_mu]/k_sigma[index_z])) * e_th

    #e_th=0.05

    #index_z = self.nbin-1
    #index_mu =self.mu_size-1
    #for index_k in range(self.k_size):
    #  print self.k[index_k,index_z,index_mu],pk_nl_th[index_k,index_z,index_mu]

    for index_z in range(2*self.nbin+1):
      for index_mu in range(self.mu_size):
        pk_nl_th[:,index_z,index_mu] *= (1.+data.mcmc_parameters['epsilon']['current']*data.mcmc_parameters['epsilon']['scale']*np.log(1. + self.k[:,index_z,index_mu]/k_sigma[index_z]) / (1. + np.log(1.+ self.k[:,index_z,index_mu]/k_sigma[index_z])) * e_th)

    #index_z = self.nbin-1
    #index_mu =self.mu_size-1
    #print ''
    #for index_k in range(self.k_size):
    #  print self.k[index_k,index_z,index_mu],pk_nl_th[index_k,index_z,index_mu]

    # Compute the beta function for nl, 
    # beta(k,z) = 1/2b(z) * d log(P_nl_th (k,z))/d log a
    #   	= -1/2b(z) *(1+z) d log(P_nl_th (k,z))/dz 
    beta_th = np.zeros((self.k_size,self.nbin,self.mu_size),'float64')
    for index_k in range(self.k_size):
      for index_z in range(self.nbin):
	  beta_th[index_k,index_z,:] = -1./(2.*self.b[index_z]) * (1.+self.z_mean[index_z]) * np.log(pk_nl_th[index_k,2*index_z+2,:]/pk_nl_th[index_k,2*index_z,:])/(self.dz)

    # Approximate the beta function without the error, and create the \tilde P_th_correction 
    #self.tilde_P_th_corr = np.zeros( (self.k_size,self.nbin,self.mu_size), 'float64')
    #for index_k in range(self.k_size):
    #  for index_z in range(self.nbin):
    #    self.tilde_P_th_corr[index_k,index_z,:] = H[2*index_z+1]/(D_A[2*index_z+1]**2) * (1. + beta_th[index_k,index_z,:]*mu[:]*mu[:])**2*data.mcmc_parameters['epsilon']['current']*data.mcmc_parameters['epsilon']['scale']*E_th[index_k,2*index_z+1,:]* pk_nl_th[index_k,2*index_z+1,:] *np.exp(-self.k[index_k,2*index_z+1,:]**2*mu[:]**2*sigma_r[index_z]**2)

    # Compute \tilde P_th(k,mu,z) = H(z)/D_A(z)^2 * (1 + beta(z,k) mu^2)^2 P_nl_th (k,z) exp(-k^2 mu^2 sigma_r^2)
    # Compute \tilde P_th(k,mu,z) = H(z)/D_A(z)^2 * (1 + beta(z,k) mu^2)^2 P_nl_th (k,z) (1 + epsilon* E(k,z) ) exp(-k^2 mu^2 sigma_r^2)
    self.tilde_P_th = np.zeros( (self.k_size,self.nbin,self.mu_size), 'float64')
    for index_k in range(self.k_size):
      for index_z in range(self.nbin):
        self.tilde_P_th[index_k,index_z,:] = H[2*index_z+1]/(D_A[2*index_z+1]**2) * (1. + beta_th[index_k,index_z,:]*mu[:]*mu[:])**2* pk_nl_th[index_k,2*index_z+1,:]*np.exp(-self.k[index_k,2*index_z+1,:]**2*mu[:]**2*sigma_r[index_z]**2)

    # Shot noise spectrum, deduced from the nuisance parameter P_shot
    self.P_shot = np.zeros( (self.nbin),'float64')
    for index_z in range(self.nbin):
      self.P_shot[index_z] = self.H_fid[2*index_z+1]/(self.D_A_fid[2*index_z+1]**2*self.b[index_z]**2)*(data.mcmc_parameters['P_shot']['current']*data.mcmc_parameters['P_shot']['scale'] + 4.*pi*r[2*index_z+1]**2*(r[2*index_z+2]-r[2*index_z])/self.n_g[index_z])
      
    #for index_z in range(self.nbin):
    #  for index_k in range(self.k_size):
    #    for index_mu in range(self.mu_size):
    #      self.tilde_P_fid[index_k,index_z,index_mu] = self.tilde_P_th[index_k,index_z,index_mu]*(1.+2.*self.alpha[index_k,2*index_z+1,index_mu])

    #index_mu = (self.mu_size-1)/2
    #index_z = 0
    #for index_k in range(1,self.k_size):
    #  print self.k_fid[index_k],self.tilde_P_th[index_k,index_z,index_mu],self.P_shot[index_z],(self.tilde_P_th[index_k,index_z,index_mu] + self.P_shot[index_z])*2.*pi/sqrt(self.k_fid[index_k]**3*self.V_survey[index_z]*self.nbin*log(self.kmax/self.kmin)),self.alpha[index_k,2.*index_z+1,index_mu]*self.tilde_P_th[index_k,index_z,index_mu]
    #print
    #print
    #index_z = self.nbin-1
    #for index_k in range(1,self.k_size):
    #  print self.k_fid[index_k],self.tilde_P_th[index_k,index_z,index_mu],self.P_shot[index_z],(self.tilde_P_th[index_k,index_z,index_mu] + self.P_shot[index_z])*2.*pi/sqrt(self.k_fid[index_k]**3*self.V_survey[index_z]*self.nbin*log(self.kmax/self.kmin)),self.alpha[index_k,2.*index_z+1,index_mu]*self.tilde_P_th[index_k,index_z,index_mu]
    #exit()
      
    # finally compute chi2, for each z_mean
    chi2 = 0.0

    mu_integrand_lo,mu_integrand_hi = 0.0,0.0
    k_integrand  = np.zeros(self.k_size,'float64')

    for index_z in range(self.nbin):
      k_integrand = self.integrand(index_z,0)
      mu_integrand_hi = np.sum((k_integrand[1:] + k_integrand[0:-1])*.5*(self.k_fid[1:] - self.k_fid[:-1]))
      for index_mu in range(1,self.mu_size):
	mu_integrand_lo = mu_integrand_hi
	mu_integrand_hi = 0
	k_integrand = self.integrand(index_z,index_mu)
	mu_integrand_hi = np.sum((k_integrand[1:] + k_integrand[0:-1])*.5*(self.k_fid[1:] - self.k_fid[:-1]))
	chi2 += (mu_integrand_hi + mu_integrand_lo)/2.*(mu[index_mu] - mu[index_mu-1])

    chi2 += (data.mcmc_parameters['epsilon']['current']*data.mcmc_parameters['epsilon']['scale'])**2

    #print chi2
    #exit()

    return - chi2/2.

  #def integrand(self,index_z,index_mu):
    #return self.k_fid[:]**2/(2.*pi)**2*self.V_survey[index_z]/2.*((self.tilde_P_th[:,index_z,index_mu] +self.tilde_P_th_corr[:,index_z,index_mu] - self.tilde_P_fid[:,index_z,index_mu])/(self.tilde_P_th[:,index_z,index_mu]+self.tilde_P_th_corr[:,index_z,index_mu]+self.P_shot[index_z]))**2
  def integrand(self,index_z,index_mu):
    return self.k_fid[:]**2/(2.*pi)**2*((self.tilde_P_th[:,index_z,index_mu] - self.tilde_P_fid[:,index_z,index_mu])**2/((2./self.V_survey[index_z])*(self.tilde_P_th[:,index_z,index_mu] + self.P_shot[index_z])**2 + (self.alpha[:,2.*index_z+1,index_mu]*self.tilde_P_th[:,index_z,index_mu])**2*self.k_fid[:]**3/2./pi**2*self.nbin*log(self.kmax/self.kmin)))

